#include <omp.h>
#include <iomanip>
#include <algorithm>

#include "inner.hpp"
#include "chromosomes.hpp"

using namespace std ;

InnerKmerExtractor::InnerKmerExtractor(int threads): threads(threads) {}

void InnerKmerExtractor::run() {
    preprocess_inner_kmers() ;
}

void dump_inner_kmers(unordered_map<uint64_t, Kmer>& kmers, Track& track) {
    nlohmann::json payload ;
    auto c = Configuration::getInstance() ;
    string p = c->workdir + "/CppInnerKmers/" + track.get_name() + ".json" ;
    std::ofstream o(p) ;
    o << "{\n" ;
    uint64_t i = 0 ;
    for (auto it = kmers.begin(); it != kmers.end(); it++) {
        if (i != 0) {
            o << ",\n" ;
        }
        o << "\"" << decode_kmer(it->first) << "\":" ;
        o << nlohmann::json(it->second).dump(4) ;
        i++ ;
    }
    o << "\n}\n" ;
}

void InnerKmerExtractor::preprocess_inner_kmers() {
    cout << "--------------------------------------------------------- " << endl ;
    cout << "Extracting inner kmers.." << endl ;
    std::unordered_map<Track, std::unordered_map<uint64_t, Kmer>> kmers ;
    for (auto it = bed_tracks.begin(); it != bed_tracks.end(); it++) {
        kmers[*it] = unordered_map<uint64_t, Kmer>() ;
    }
    int n = 0 ;
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < bed_tracks.size(); i++) {
        if (bed_tracks[i].svtype == SVTYPE_DEL) {
            kmers[bed_tracks[i]] = extract_deletion_inner_kmers(bed_tracks[i]) ;
        }
        if (bed_tracks[i].svtype == SVTYPE_INS) {
            kmers[bed_tracks[i]] = extract_insertion_inner_kmers(bed_tracks[i]) ;
        }
        n += 1 ;
        dump_inner_kmers(kmers[bed_tracks[i]], bed_tracks[i]) ;
        if (n % 100 == 0) {
            cout << "Loaded " << std::setw(6) << std::left << n << " tracks approximately.." << '\r' << flush ;
        }
    }
    cout << endl ;
    cout << "Merging inner kmers.." << endl ;
    int i = 0 ;
    auto it = kmers.begin() ;
    while (it != kmers.end()) {
        for (auto kmer = it->second.begin(); kmer != it->second.end(); kmer++) {
            if (inner_kmers.find(kmer->first) == inner_kmers.end()) {
                inner_kmers[kmer->first] = kmer->second ;
                inner_kmers[kmer->first].loci.clear() ;
                inner_kmers[kmer->first].tracks.clear() ;
            }
            for (auto track = kmer->second.tracks.begin(); track != kmer->second.tracks.end(); track++) {
                if (inner_kmers[kmer->first].tracks.find(track->first) == inner_kmers[kmer->first].tracks.end()) {
                    inner_kmers[kmer->first].tracks[track->first] = 0 ;
                }
                inner_kmers[kmer->first].tracks[track->first] += track->second ;
            }
            auto& loci = inner_kmers[kmer->first].loci ;
            loci.insert(loci.begin(), kmer->second.loci.begin(), kmer->second.loci.end()) ;
        }
        if (i % 100 == 0) {
            cout << "Loaded " << std::setw(6) << std::left << i << " tracks approximately.." << '\r' << flush ;
        }
        i += 1 ;
        it = kmers.erase(it) ;
    }
    cout << endl ;
    cout << inner_kmers.size() << " inner kmers.." << endl ;
}

unordered_map<uint64_t, Kmer> InnerKmerExtractor::extract_deletion_inner_kmers(Track& track) {
    SimpleTrack _track ;
    _track.end = track.end ;
    _track.begin = track.begin ;
    _track.chrom = track.chrom ;
    _track.svtype = track.svtype ;
    unordered_map<uint64_t, Kmer> kmers ;
    int l = track.svlen + 32 + 32 ;
    char* seq = (char*) malloc(sizeof(char) * (l + 1)) ;
    strncpy(seq, chromosome_seqs[get_chromosome_name(track.chrom)] + track.begin - 32, l) ;
    seq[l] = '\0' ;
    KmerIterator it(seq, 0, l - 32 + 1, track.begin - 32, ITER_MODE_READ) ;
    while (it) {
        auto canon = encode_kmer(canonicalize(decode_kmer(it->kmer))) ;
        if (kmers.find(canon) == kmers.end()) {
            kmers[canon] = Kmer(canon, KMER_TYPE_INNER) ;
            kmers[canon].trend = false ;
        }
        // loci will be found by scanner
        auto& tracks = kmers[canon].tracks ;
        if (tracks.find(_track) == tracks.end()) {
            tracks[_track] = 0 ;
        }
        tracks[_track] += 1 ;
        it++ ;
    }
    free(seq) ;
    //cout << inner_kmers.size() << " inner kmers." << endl ;
    return kmers ;
}

unordered_map<uint64_t, Kmer> InnerKmerExtractor::extract_insertion_inner_kmers(Track& track) {
    SimpleTrack _track ;
    _track.end = track.end ;
    _track.begin = track.begin ;
    _track.chrom = track.chrom ;
    _track.svtype = track.svtype ;
    int offset = 250 ;
    int l = 64 ;
    unordered_map<uint64_t, Kmer> breakpoint_kmers ;
    char* breakpoint_seq = (char*) malloc(sizeof(char) * (l + 1)) ;
    strncpy(breakpoint_seq, chromosome_seqs[get_chromosome_name(track.chrom)] + track.begin - 32, 64) ;
    breakpoint_seq[l] = '\0' ;
    KmerIterator it(breakpoint_seq, 0, 32 + 1, track.begin - 32, ITER_MODE_READ) ; 
    while (it) {
        auto canon = encode_kmer(canonicalize(decode_kmer(it->kmer))) ;
        if (breakpoint_kmers.find(canon) == breakpoint_kmers.end()) {
            breakpoint_kmers[canon] = Kmer(canon, KMER_TYPE_INNER) ;
            breakpoint_kmers[canon].trend = false ;
        }
        // loci will be found by scanner
        auto& tracks = breakpoint_kmers[canon].tracks ;
        if (tracks.find(_track) == tracks.end()) {
            tracks[_track] = 0 ;
        }
        tracks[_track] += 1 ;
        it++ ;
    }
    free(breakpoint_seq) ;
    //cout << breakpoint_kmers.size() << " breakpoint kmers." << endl ;
    // inner sequence
    l = track.seq.length() + offset + offset ;
    char* inner_seq = (char*) malloc(sizeof(char) * (l + 1)) ;
    strncpy(inner_seq, chromosome_seqs[get_chromosome_name(track.chrom)] + track.begin - offset, offset) ;
    strncpy(inner_seq + offset, track.seq.c_str(), track.seq.length()) ;
    strncpy(inner_seq + offset + track.seq.length(), chromosome_seqs[get_chromosome_name(track.chrom)] + track.end, offset) ;
    inner_seq[l] = '\0' ;
    for (int i = 0; i < strlen(inner_seq); i++) {
        inner_seq[i] = toupper(inner_seq[i]) ;
    }
    unordered_map<uint64_t, Kmer> inner_kmers ;
    it = KmerIterator(inner_seq, offset + 32, l - offset - 64 + 1, track.begin - 32, ITER_MODE_READ) ;
    while (it) {
        auto canon = encode_kmer(canonicalize(decode_kmer(it->kmer))) ;
        if (breakpoint_kmers.find(canon) != breakpoint_kmers.end()) {
            breakpoint_kmers.erase(breakpoint_kmers.find(canon)) ;
            it++ ;
            continue ;
        }
        if (inner_kmers.find(canon) == inner_kmers.end()) {
            inner_kmers[canon] = Kmer(canon, KMER_TYPE_INNER) ;
            inner_kmers[canon].trend = true ;
        }
        //cout << decode_kmer(it->left) << endl ;
        //cout << decode_kmer(it->kmer) << endl ;
        //cout << decode_kmer(it->right) << endl ;
        if (!it.has_left || !it.has_right) {
            cout << track.get_name() << endl ;
            cout << it.seq << endl ;
            cout << it.i << endl ;
        }
        assert(it.has_left && it.has_right) ;
        Locus locus({track.chrom, track.begin, LOCUS_TYPE_INNER, it->left, it->right, it->gc / 5}) ;
        inner_kmers[canon].loci.push_back(locus) ;
        auto& tracks = inner_kmers[canon].tracks ;
        if (tracks.find(_track) == tracks.end()) {
            tracks[_track] = 0 ;
        }
        tracks[_track] += 1 ;
        it++ ;
    }
    free(inner_seq) ;
    inner_kmers.insert(breakpoint_kmers.begin(), breakpoint_kmers.end()) ;
    //cout << inner_kmers.size() << " inner kmers." << endl ;
    return inner_kmers ;
}
