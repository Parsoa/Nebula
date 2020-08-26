#include <omp.h>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "logger.hpp"
#include "junction.hpp"

using namespace std ;

JunctionKmerExtractor::JunctionKmerExtractor(int threads): threads(threads) {}

void JunctionKmerExtractor::run() {
    cout << "--------------------------------------------------------- " << endl ;
    load_tracks() ;
    preprocess_junction_kmers() ;
}

void JunctionKmerExtractor::load_tracks() {
    auto c = Configuration::getInstance() ;
    // checks if BAM files exist
    for (auto bam_file = c->bam.begin(); bam_file != c->bam.end(); bam_file++) {
        samFile *in = sam_open((*bam_file).c_str(), "r") ;
        if (in == NULL) {
            cerr << "BAM file " << *bam_file << " does not exist. Aborting.." << endl ;
            exit(0) ;
        }
    }
    int i = 0 ;
    for (auto it = c->bed.begin(); it != c->bed.end(); it++) {
        files.push_back(std::vector<samFile*>(threads)) ;
        tracks.push_back(load_tracks_from_file_as_dict(*it)) ;
        headers.push_back(std::vector<bam_hdr_t*>(threads)) ;
        indices.push_back(std::vector<hts_idx_t*>(threads)) ;
        iterators.push_back(std::vector<hts_itr_t*>(threads)) ;
        for (int j = 0; j < threads; j++) {
            files[i][j] = sam_open(c->bam[i].c_str(), "r") ;
            headers[i][j] = sam_hdr_read(files[i][j]) ;
            indices[i][j] = sam_index_load(files[i][j], c->bam[i].c_str()) ;
        }
        i++ ;
    }
}

void dump_junction_kmers(unordered_map<uint64_t, Kmer>& kmers, Track& track) {
    nlohmann::json payload ;
    auto c = Configuration::getInstance() ;
    string p = c->workdir + "/CppJunctionKmers/" + track.get_name() + ".json" ;
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

void JunctionKmerExtractor::preprocess_junction_kmers() {
    cout << "Extracting junction kmers.." << endl ;
    std::unordered_map<Track, std::unordered_map<uint64_t, Kmer>> kmers ;
    for (auto it = bed_tracks.begin(); it != bed_tracks.end(); it++) {
        kmers[*it] = unordered_map<uint64_t, Kmer>() ;
    }
    int n = 0 ;
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < bed_tracks.size(); i++) {
        auto track = bed_tracks[i] ;
        //if (track.get_name() != "INS@chr21_25310099_25310100") {
        //    continue ;
        //}
        int t = omp_get_thread_num() ;
        //cout << bed_tracks[i].get_name() << std::endl ;
        std::unordered_map<uint64_t, Kmer> absent_kmers ;
        std::unordered_map<uint64_t, Kmer> present_kmers ;
        bool first_batch = true ;
        for (int j = 0; j < tracks.size(); j++){
            auto _kmers = extract_junction_kmers(track, t, j) ;
            //cout << _kmers.size() << endl ;
            if (tracks[j].find(track) == tracks[j].end() || tracks[j].find(track)->first.genotype == "0/0" || tracks[j].find(track)->first.genotype == "./.") {
                //cout << "absent" << endl ;
                absent_kmers.insert(_kmers.begin(), _kmers.end()) ;
            } else {
                //cout << "present" << endl ;
                if (first_batch) {
                    first_batch = false ;
                    present_kmers.insert(_kmers.begin(), _kmers.end()) ;
                } else {
                    auto kmer = present_kmers.begin() ;
                    while (kmer != present_kmers.end()) {
                        if (_kmers.find(kmer->first) == _kmers.end()) {
                            kmer = present_kmers.erase(kmer) ;
                        } else {
                            kmer++ ;
                        }
                    }
                } 
            }
        }
        //cout << "Opposing kmers: " << absent_kmers.size() << endl ;
        //cout << "Supporting kmers:" << present_kmers.size() << endl ;
        for (auto kmer = present_kmers.begin(); kmer != present_kmers.end(); kmer++){
            if (absent_kmers.find(kmer->first) == absent_kmers.end()) {
                kmers[track][kmer->first] = kmer->second ;
            }
        }
        dump_junction_kmers(kmers[track], track) ;
        n += 1 ;
        if (n % 100 == 0) {
            cout << "Loaded " << std::setw(6) << std::left << n << " tracks approximately.." << '\r' << flush ;
        }
    }
    cout << endl ;
    cout << "Merging junction kmers.." << endl ;
    for (auto track = kmers.begin(); track != kmers.end(); track++) {
        for (auto kmer = track->second.begin(); kmer != track->second.end(); kmer++) {
            if (junction_kmers.find(kmer->first) == junction_kmers.end()) {
                junction_kmers[kmer->first] = kmer->second ;
                junction_kmers[kmer->first].loci.clear() ;
                junction_kmers[kmer->first].tracks.clear() ;
            }
            for (auto track = kmer->second.tracks.begin(); track != kmer->second.tracks.end(); track++) {
                if (junction_kmers[kmer->first].tracks.find(track->first) == junction_kmers[kmer->first].tracks.end()) {
                    junction_kmers[kmer->first].tracks[track->first] = 0 ;
                }
                junction_kmers[kmer->first].tracks[track->first] += track->second ;
            }
            auto& loci = junction_kmers[kmer->first].loci ;
            loci.insert(loci.begin(), kmer->second.loci.begin(), kmer->second.loci.end()) ;
        }
    }
    // Check kmer consistency
    int e = 0 ;
    auto kmer = junction_kmers.begin() ;
    while (kmer != junction_kmers.end()) {
        int s = 0 ; 
        for (auto __track = kmer->second.tracks.begin(); __track != kmer->second.tracks.end(); __track++) {
            s += __track->second ;
        }
        // This can happen when kmer is very common, avoid synchornization and ignore kmer
        if (s != kmer->second.loci.size()) {
            e += 1 ;
            kmer = junction_kmers.erase(junction_kmers.find(kmer->first)) ;
            continue ;
        }
        kmer++ ;
    }
    cout << "Erased " << e << " corrupted kmers." << endl ;
    cout << junction_kmers.size() << " junction kmers.." << endl ;
}

unordered_map<uint64_t, Kmer> JunctionKmerExtractor::extract_junction_kmers(Track track, int t, int j) {
    SimpleTrack _track ;
    _track.end = track.end ;
    _track.begin = track.begin ;
    _track.chrom = track.chrom ;
    _track.svtype = track.svtype ;
    //
    hts_itr_t* iter = iterators[j][t] ;
    int id = bam_name2id(headers[j][t], get_chromosome_name(track.chrom).c_str()) ;
    iter = sam_itr_queryi(indices[j][t], id, track.begin - 100, track.end + 100) ;
    bam1_t* read = bam_init1() ;
    unordered_map<uint64_t, Kmer> kmers ;
    int m = 0;
    int n = 0 ;
    while (sam_itr_next(files[j][t], iter, read) >= 0) {
        if (read->core.pos < track.begin - 100 || read->core.pos > track.end + 100){
            continue ;
        }
        n += 1 ;
        int offset = 0 ;
        uint32_t* cigar = bam_get_cigar(read) ;
        std::vector<std::pair<int, int>> clips ;
        bool left_clipped = false ;
        bool right_clipped = false ;
        for (int i = 0; i < read->core.n_cigar; i++) {
            int l = bam_cigar_oplen(cigar[i]) ;
            if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
                clips.push_back(std::make_pair(offset, offset + l)) ;
                if (i == 0) {
                    left_clipped = true ;
                } else {
                    right_clipped = true ;
                }
            }
            else if (bam_cigar_op(cigar[i]) == BAM_CINS || bam_cigar_op(cigar[i]) == BAM_CDEL) {
                if (l >= 0.9 * double(track.svlen) && l <= 1.1 * double(track.svlen)) {
                    clips.push_back(std::make_pair(offset, offset + l)) ;
                }
            }
            offset += l ;
        }
        //for (auto p = clips.begin(); p != clips.end(); p++) {
        //    cout << p->first << "," << p->second << " " ;
        //}
        //cout << endl ;
        if (right_clipped && left_clipped) {
            continue ;
        }
        uint32_t l = read->core.l_qseq ; //length of the read
        char* seq = (char*) malloc(l + 1) ;
        uint8_t *q = bam_get_seq(read) ; //quality string
        for (int i = 0; i < l; i++){
            seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        seq[l] = '\0' ;
        if (strlen(seq) < 40) {
            continue ;
        }
        KmerIterator it(seq, 0, l - 32 + 1, 0, ITER_MODE_READ) ;
        while (it) {
            if (it->kmer != 0 && is_clipped(clips, std::make_pair(it.i - 1, it.i + 31))) {
                auto canon = encode_kmer(canonicalize(decode_kmer(it->kmer))) ;
                m += 1 ;
                //cout << it.i - 1 << " " << it.i + 31 << endl ;
                if (kmers.find(canon) == kmers.end()) {
                    kmers[canon] = Kmer(canon, KMER_TYPE_JUNCTION) ;
                    kmers[canon].seq = canon ;
                    kmers[canon].trend = true ;
                }
                // AH AH, RACE Condition
                // only one junction loci per track
                // TODO: what if multiple tracks insert at the same size?
                if (kmers[canon].loci.size() == 0) {
                    Locus locus {track.chrom, track.begin, LOCUS_TYPE_JUNCTION, it->left, it->right, it->gc} ;
                    auto d = kmers[canon] ;
                    kmers[canon].loci.push_back(locus) ;
                } else {
                    Locus& locus = kmers[canon].loci[0] ;
                    if (locus.left == 0) {
                        locus.left = it->left ;
                    }
                    if (locus.right == 0) {
                        locus.right = it->right ;
                    }
                }
                if (kmers[canon].tracks.find(_track) == kmers[canon].tracks.end()) {
                    kmers[canon].tracks[_track] = 1 ;
                }
            }
            it++ ;
        }
    }
    //cout << m << endl ;
    //cout << n << " reads " << endl ;
    //cout << kmers.size() << endl ;
    bam_destroy1(read) ;
    return kmers ;
}

int JunctionKmerExtractor::overlap(std::pair<int, int> a, std::pair<int, int> b) {
    return std::max(0, std::min(a.second, b.second) - std::max(a.first, b.first)) ;
}

bool JunctionKmerExtractor::is_clipped(std::vector<std::pair<int, int>> clips, std::pair<int, int> kmer) {
    for (auto clip = clips.begin(); clip != clips.end(); clip++) {
        if (overlap(*clip, kmer) >= 10) {
            return true ;
        }
    }
    return false ;
}

