#include <omp.h>
#include <ctime>
#include <string>
#include <thread>
#include <locale>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream> 
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <pthread.h>
#include <unordered_map>

#include "json.hpp"
#include "bed_utils.hpp"
#include "kmer_utils.hpp"

using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

int KMER_TYPE ;
std::unordered_map<uint64_t, std::vector<uint64_t>> counts ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

#define LOCUS_TYPE_JUNCTION uint8_t(0)
#define LOCUS_TYPE_INNER    uint8_t(1)
#define LOCUS_TYPE_REF      uint8_t(2)

struct Locus {
    uint8_t chrom ;
    uint32_t position ;
    uint8_t type ; 
    uint64_t left ;
    uint64_t right ;
    uint16_t gc ;
} ;

struct Kmer {
    uint32_t count ;
    uint32_t total ;
    uint32_t reference ;
    std::unordered_map<string, uint8_t> tracks ;
    std::vector<Locus> loci ;
    std::vector<Locus> filtered_loci ;
    std::vector<Locus> junction_loci ;
    bool inverse ;
    Kmer(): inverse(false), count(0), total(0), reference(0) {}
} ;

void parse_locus_name(string name, Locus &locus) {
    stringstream ss(name) ;
    vector<string> tokens ;
    string token ;
    while (getline(ss, token, '_')) {
        tokens.push_back(token) ;
    }
    locus.type = tokens[0] == "junction" ? LOCUS_TYPE_JUNCTION : LOCUS_TYPE_INNER ;
    locus.chrom = get_chromosome_index(tokens[1].substr(tokens[1].find("@") + 1)) ;
    locus.position = std::stoi(tokens[2]) ;
}

string get_locus_name(Locus locus) {
    if (locus.type == LOCUS_TYPE_INNER) {
        return "inside_" + get_chromosome_name(locus.chrom) + "_" + std::to_string(locus.position) ;
    }
    if (locus.type == LOCUS_TYPE_JUNCTION) {
        return "junction_" + get_chromosome_name(locus.chrom) + "_" + std::to_string(locus.position) ;
    }
    return get_chromosome_name(locus.chrom) + "_" + std::to_string(locus.position) ;
}

std::unordered_map<string, Locus> cast_loci_map(vector<Locus>& _loci) {
    unordered_map<string, Locus> loci ;
    for (auto locus = _loci.begin(); locus != _loci.end(); locus++) {
        loci[get_locus_name(*locus)] = *locus ;
    }
    return loci ;
}

unordered_map<string, Locus> jsonify_loci(vector<Locus> loci) {
    unordered_map<string, Locus> _loci ;
    for (auto locus = loci.begin(); locus != loci.end(); locus++) {
        _loci[get_locus_name(*locus)] = *locus ;
    }
    return _loci ;
}

void to_json(nlohmann::json& j, const Kmer& k) {
    j = nlohmann::json{{"count", k.count}, {"reference", k.reference}, {"loci", jsonify_loci(k.loci)}, {"junction_loci", jsonify_loci(k.junction_loci)}, {"filtered_loci", jsonify_loci(k.filtered_loci)}, {"tracks", k.tracks}, {"inverse", k.inverse}} ;
}

void to_json(nlohmann::json& j, const Locus& l) {
    nlohmann::json masks({}) ;
    if (l.left != 0 ) {
        masks[decode_kmer(l.left)] = 0 ;
    }
    if (l.right != 0) {
        masks[decode_kmer(l.right)] = 0 ;
    }
    j = nlohmann::json{{"chrom", get_chromosome_name(l.chrom)}, {"position", l.position}, {"masks", masks}, {"type", l.type}} ;
}

void from_json(const nlohmann::json& j, Locus& l) {
    nlohmann::json masks = j.at("masks") ;
    if (masks.find("left") != masks.end()) {
        l.left = encode_kmer(masks.at("left")) ;
    } else {
        l.left = 0 ;
    }
    if (masks.find("right") != masks.end()) {
        l.right = encode_kmer(masks.at("right")) ;
    } else {
        l.right = 0 ;
    }
}

void from_json(const nlohmann::json& j, Kmer& k) {
    for (auto locus = j.at("loci").begin(); locus != j.at("loci").end(); locus++) {
        Locus l = locus.value().get<Locus>() ;
        parse_locus_name(locus.key(), l) ;
        k.loci.push_back(l) ;
    }
    if (j.find("reference") != j.end()) {
        k.reference = j["reference"] ;
    }
}

std::unordered_map<uint64_t, Kmer> kmers ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

#define MASK 6
void process_read(string chrom, const char* seq, unordered_map<uint64_t, Kmer>& _kmers) {
    int gc = 0 ;
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    uint64_t l = strlen(seq) ;
    //cout << "Scanning " << chrom << ".." << endl ;
    for (uint64_t i = 250 - 16; i < l - 500; i++) {
        if (i == 250 - 16) {
            left = encode_kmer(seq + i - 32) ;
            k = encode_kmer(seq + i) ;
            right = encode_kmer(seq + i + 32) ;
            gc = calc_gc_content(seq + i) ;
        } else {
            k = k << 2 ;
            k += (seq[i + 31] & MASK) >> 1 ;
            right = right << 2 ;
            right += (seq[i + 32 + 31] & MASK) >> 1 ;
            left = left << 2 ;
            left += (seq[i - 1] & MASK) >> 1 ;
            if (seq[i - (250 - 16)] == 'C' || seq[i - (250 - 16)] == 'G') {
                gc -= 1 ;
            }
            if (seq[i + (250 + 32)] == 'C' || seq[i + (250 + 32)] == 'G') {
                gc += 1 ;
            }
        }
        if (kmers.find(k) != kmers.end()) {
            if (_kmers.find(k) == _kmers.end()) {
                _kmers[k] = Kmer() ;
                _kmers[k].count = 0 ;
            }
            string name = chrom + "@" + to_string(i) ;
            Locus locus({get_chromosome_index(chrom), uint32_t(i), LOCUS_TYPE_REF, left, right, gc}) ;
            _kmers[k].loci.push_back(locus) ;
            _kmers[k].count ++ ;
        } else {
            uint64_t rc = encoded_reverse_complement(k) ;
            if (kmers.find(rc) != kmers.end()) {
                if (_kmers.find(rc) == _kmers.end()) {
                    _kmers[rc] = Kmer() ;
                    _kmers[rc].count = 0 ;
                }
                string name = chrom + "@" + to_string(i) ;
                Locus locus({get_chromosome_index(chrom), uint32_t(i), LOCUS_TYPE_REF, left, right, gc}) ;
                _kmers[rc].loci.push_back(locus) ;
                _kmers[rc].count ++ ;
            }
        }
        if (i % 1000000 == 0) {
            cout << chrom << " progress " << float(i) / float(l) << ".." << endl ;
        }
    }
}

// ============================================================================= \\
// ================================ FASTA Files ================================ \\
// ============================================================================= \\

vector<string> chromosomes ;
unordered_map<string, char*> chromosome_seqs ;

int get_reference_size(ifstream &fasta_file) {
    fasta_file.seekg(0, ios_base::end) ;
    int l = fasta_file.tellg() ;
    fasta_file.seekg(0, ios_base::beg) ;
    return l ;
}

void load_chromosomes(string path) {
    // assume human genome length
    ifstream fasta_file ;
    fasta_file.open(path, ios::binary) ;
    int l = get_reference_size(fasta_file) ;
    // maximum size of a chromosome, kinda arbitrary
    char* buffer = (char*) malloc(sizeof(char) * 300000000) ;
    // read all of file
    int state ;
    uint64_t n = 0 ;
    std::string line ;
    std::getline(fasta_file, line) ;
    while (true) {
        if (line.substr(0, 4) == ">chr") {
            int l = line.length() ;
            if (l == 6 || l == 5) {
                if (line[4] == 'X' || line[4] == 'Y' || (line[4] >= '1' && line[4] <= '9')) {
                    string chrom = line.substr(1, l - 1) ;
                    cout << "Collecting " << chrom << ".." << endl ;
                    while(std::getline(fasta_file, line)) {
                         if (line[0] == '>') {
                             break ;
                         }
                         for (int i = 0; i < line.length(); i++) {
                            line[i] = toupper(line[i]) ;
                         }
                         memcpy(buffer + n, line.c_str(), line.length()) ;
                         n += line.length() ;
                    }
                    buffer[n] = '\0' ;
                    cout << "Extracted " << chrom << " with " << n << " bases." << endl ;
                    char* s = (char*) malloc(sizeof(char) * (n + 1)) ;
                    memcpy(s, buffer, n + 1) ;
                    chromosomes.push_back(chrom) ;
                    chromosome_seqs[chrom] = s ;
                    n = 0 ;
                    continue ;
                }
            }
        } 
        if (!std::getline(fasta_file, line)) {
            break ;
        }
    }
    free(buffer) ;
}


unordered_map<uint64_t, Kmer> scan_chromosome(string chrom) {
    unordered_map<uint64_t, Kmer> _kmers ;
    process_read(chrom, chromosome_seqs[chrom], _kmers) ;
    free(chromosome_seqs[chrom]) ;
    return _kmers ;
}

void scan_reference(int threads) {
    int m = 0 ;
    int n = chromosomes.size() ;
    vector<unordered_map<uint64_t, Kmer>> batches(threads) ; 
    while (m < n) {
        int p = m ;
        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < threads; i++) {
            if (m + i < n) {
                batches[i] = scan_chromosome(chromosomes[m + i]) ;
            }
        }
        m += threads ;
        if (m >= chromosomes.size()) {
            m = n ;
        }
        // merge batches
        for (int i = p; i < m; i++) {
            cout << "Merging " << chromosomes[i] << ", " << batches[i - p].size() << " matches.." << endl ;
            for(auto it = batches[i - p].begin(); it != batches[i - p].end(); it++) {
                kmers[it->first].count += it->second.count ;
                kmers[it->first].loci.insert(kmers[it->first].loci.end(), it->second.loci.begin(), it->second.loci.end()) ;
            }
            batches[i - p].clear() ;
        }
    }
    cout << "Reference scan completed." << endl ;
}

// ============================================================================= \\
// ================================= Filtering ================================= \\
// ============================================================================= \\

vector<uint64_t> find_interest_masks(Kmer& k) {
    std::vector<uint64_t> interset_masks ;
    for (auto locus = k.loci.begin(); locus != k.loci.end(); locus++) {
        if (locus->type != LOCUS_TYPE_REF) {
            if (locus->left != 0) {
                interset_masks.push_back(locus->left) ;
            }
            if (locus->right != 0) {
                interset_masks.push_back(locus->right) ;
            }
        }
        // ref loci that fall inside event
        else {
            for (auto track = k.tracks.begin(); track != k.tracks.end(); track++) {
                Track t = parse_track_name(track->first) ;
                if (locus->chrom == t.chrom) {
                    if (locus->position >= t.begin - 32 && locus->position < t.end) {
                        if (locus->left != 0) {
                            interset_masks.push_back(locus->left) ;
                        }
                        if (locus->right != 0) {
                            interset_masks.push_back(locus->right) ;
                        }
                    }
                } 
            }
        }
    }
    return interset_masks ;
}

//TODO: does this even help?
// removes kmers that have loci too close the breakpoints
bool is_kmer_returning(Kmer& kmer) {
    for(auto track = kmer.tracks.begin(); track != kmer.tracks.end(); track++) {
        Track t = parse_track_name(track->first) ;
        for(auto locus = kmer.loci.begin(); locus != kmer.loci.end(); locus++) {
            // inner kmers for deletion only have ref loci
            if (locus->type == LOCUS_TYPE_REF && !(KMER_TYPE == LOCUS_TYPE_INNER && t.svtype == SVTYPE_DEL)) {
                if (abs(int(locus->position) - int(t.end)) < 64) {
                    return true ;
                }
                if (abs(int(locus->position) - int(t.begin)) < 64) {
                    return true ;
                }
            }
        }
    }
    return false ;
}

void filter_kmers() {
    cout << "Filtering " << kmers.size() << " kmers based on ref count.." << endl ;
    auto kmer = kmers.begin() ;
    while (kmer != kmers.end()) {
        Kmer& k = kmer->second ;
        k.reference = k.count ;
        if (k.count > 3) {
            kmer = kmers.erase(kmer) ;
        } else {
            kmer++ ;
        }
    }
    //TODO: Parallelize this
    cout << "Filtering " << kmers.size() << " kmers based on masks.." << endl ;
    kmer = kmers.begin() ;
    while (kmer != kmers.end()) {
        Kmer& k = kmer->second ;
        std::vector<Locus> loci(k.loci) ;
        std::vector<uint64_t> interset_masks = find_interest_masks(k) ;
        // number of loci before filtering
        int l_1 = k.loci.size() ;
        auto locus = k.loci.begin() ;
        while (locus != k.loci.end()) {
            bool found = false ;
            for (auto m = interset_masks.begin(); m != interset_masks.end(); m++) {
                if (locus->left != 0) {
                    if (is_canonical_subsequence(locus->left, *m)) {
                        found = true ;
                        break ;
                    }
                }
                if (locus->right != 0) {
                    if (is_canonical_subsequence(locus->right, *m)) {
                        found = true ;
                        break ;
                    }
                }
            }
            // doesn't have any shared masks, filter
            if (not found) {
                k.filtered_loci.push_back(*locus) ;
                locus = k.loci.erase(locus) ;
                continue ;
            }
            locus++ ;
        }
        // number of loci after filtering
        int l_2 = k.loci.size() ;
        // loci with less than two masks exist
        auto a = find_if(k.loci.begin(), k.loci.end(), [](Locus l) {
            return l.left == 0 || l.right == 0 ;
        }) ;
        //non-junction loci exists
        auto b = find_if(k.loci.begin(), k.loci.end(), [](Locus l) {
            return l.type == LOCUS_TYPE_REF ;
        }) ;
        // won't happen for inner kmers
        if (a != k.loci.end()) { // junction loci exist wtih less than two masks
            if (l_1 != l_2) { // some ref loci were filtered
                if (b == k.loci.end()) { // all ref loci were filtered
                    // because ref loci will have masks, count them instead and subtract from total. May overcount.
                    // Counting with one mask will undercount
                    k.loci = loci ;
                    auto locus = k.loci.begin() ; // count non-junction loci only
                    while (locus != k.loci.end()) {
                        if (locus->type == KMER_TYPE) {
                            k.junction_loci.push_back(*locus) ;
                            locus = k.loci.erase(locus) ;
                        } else {
                            locus++ ;
                        }
                    }
                    k.inverse = true ;
                } else { // some ref loci remain
                    kmer = kmers.erase(kmer) ;
                    continue ;
                }
            } else { // no ref loci was filtered, so we need to count every loci, ignore masks
                for (auto locus = k.loci.begin(); locus != k.loci.end(); locus++){
                    locus->left = 0 ;
                    locus->right = 0 ;
                }
            }
        }
        kmer++ ;
    }
    //remove returning kmers
    cout << "Filtering " << kmers.size() << " kmers based on coordinates.." << endl ;
    kmer = kmers.begin() ;
    while (kmer != kmers.end()) {
        Kmer& k = kmer->second ;
        if (is_kmer_returning(k)) {
            kmer = kmers.erase(kmer) ;
        } else {
            kmer++ ;
        }
    }
    cout << "Remaining " << kmers.size() << endl ;
}

void output_kmers(string path) {
    nlohmann::json payload ;
    cout << "Dumping kmer counts..." << endl ;
    string p = KMER_TYPE == LOCUS_TYPE_JUNCTION ? path + "/FilterJunctionKmersJob/kmers_cpp.json" : path + "/FilterInnerKmersJob/kmers_cpp.json" ;
    cout << p << endl ;
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

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int load_kmers(string path, int threads) {
    cout << "Loading kmers from " << path << ".." << endl ;
    ifstream json_file(KMER_TYPE == LOCUS_TYPE_JUNCTION ? path + "/ExtractJunctionKmersJob/batch_merge.json" : path + "/ExtractInnerKmersJob/batch_merge.json") ;
    nlohmann::json index_json ;
    json_file >> index_json ;
    std::vector<std::unordered_map<uint64_t, Kmer>> _kmers(threads) ;
    int n = 0 ;
    std::vector<string> filenames ;
    std::vector<string> tracknames ;
    for (nlohmann::json::iterator it = index_json.begin(); it != index_json.end(); ++it) {
        filenames.push_back(it.value()) ;
        tracknames.push_back(it.key()) ;
    }
    int l = filenames.size() ;
    ofstream num_file(path + "/ExtractInnerKmersJob/num_cpp.txt") ;
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < l; i++) {
        // load kmers:
        ifstream track_file(filenames[i]) ;
        nlohmann::json track_json ;
        track_file >> track_json ;
        // tmp, for inner kmers only
        int t = omp_get_thread_num() ;
        std::unordered_map<uint64_t, Kmer> __kmers(threads) ;
        for (nlohmann::json::iterator kmer = track_json.begin(); kmer != track_json.end(); kmer++) {
            uint64_t k = encode_kmer(canonicalize(kmer.key()).c_str()) ;
            if (KMER_TYPE == LOCUS_TYPE_INNER) {
                assert(__kmers.find(k) == __kmers.end()) ;
                __kmers[k] = kmer.value().get<Kmer>() ;
                __kmers[k].count = 0 ;
                __kmers[k].tracks[tracknames[i]] = 1 ;
            } else {
                _kmers[t][k] = kmer.value().get<Kmer>() ;
                _kmers[t][k].count = 0 ;
                _kmers[t][k].tracks[tracknames[i]] = 1 ;
            }
        }
        if (KMER_TYPE == LOCUS_TYPE_INNER) {
            assert(__kmers.size() == track_json.size()) ;
        }
        if (KMER_TYPE == LOCUS_TYPE_INNER) {
            std::vector<std::unordered_map<uint64_t, Kmer>> _ref(6) ;
            int m = 0 ;
            for (auto kmer = __kmers.begin(); kmer != __kmers.end(); kmer++) {
                if (kmer->second.reference <= 5) {
                    _ref[kmer->second.reference][kmer->first] = kmer->second ;
                    m += 1 ;
                } else {
                    num_file << tracknames[i] << " skipping " << decode_kmer(kmer->first) << endl ;
                }
            }
            int n = 0 ;
            for (int j = 0; j <= 5; j++) {
                for (auto kmer = _ref[j].begin(); kmer != _ref[j].end(); kmer++) {
                    if (n == 50) {
                        break ;
                    } else {
                        if (_kmers[t].find(kmer->first) == _kmers[t].end()) {
                            _kmers[t][kmer->first] = kmer->second ;
                        } else {
                            _kmers[t][kmer->first].loci.insert(_kmers[t][kmer->first].loci.end(), kmer->second.loci.begin(), kmer->second.loci.end()) ;
                            _kmers[t][kmer->first].tracks.insert(kmer->second.tracks.begin(), kmer->second.tracks.end()) ;
                        }
                        n += 1 ;
                    }
                }
            }
        }
        track_file.close() ;
        n += 1 ;
        if (t == 0) {
            if (n % 100 == 0) {
                cout << "Loaded " << n << " tracks approximately.." << '\r' << flush ;
            }
        }
    }
    cout << endl ;
    json_file.close() ;
    for (int i = 0; i < threads; i++) {
        for (auto kmer = _kmers[i].begin(); kmer != _kmers[i].end(); kmer++) {
            if (kmers.find(kmer->first) == kmers.end()) {
                kmers[kmer->first] = kmer->second ;
            } else {
                kmers[kmer->first].loci.insert(kmers[kmer->first].loci.end(), kmer->second.loci.begin(), kmer->second.loci.end()) ;
                kmers[kmer->first].tracks.insert(kmer->second.tracks.begin(), kmer->second.tracks.end()) ;
            }
        }
    }
    cout << "Loaded " << kmers.size() << " kmers." << endl ;
    return 0 ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

unordered_map<uint64_t, Kmer> extract_insertion_kmers(Track& track) {
    unordered_map<uint64_t, Kmer> kmers ;
    if (track.svlen < 96) {
        return kmers ;
    }
    int l = track.svlen + 250 + 250 ;
    char* seq = (char*) malloc(sizeof(char) * (l + 1)) ;
    strncpy(seq, chromosome_seqs[get_chromosome_name(track.chrom)] + track.begin - 250, 250) ;
    strcpy(seq + 250, track.seq.c_str()) ;
    strncpy(seq + 250 + track.seq.length(), chromosome_seqs[get_chromosome_name(track.chrom)] + track.end, 250) ;
    KmerIterator it(seq, 250, l - 250) ;
    return kmers ;
}

unordered_map<uint64_t, Kmer> extract_deletion_kmers(Track& track) {
    unordered_map<uint64_t, Kmer> kmers ;
    char* seq = (char*) malloc(sizeof(char) * (64 + 1)) ;
    strncpy(seq, chromosome_seqs[get_chromosome_name(track.chrom)] + track.begin - 32, 32) ;
    strncpy(seq + 32, chromosome_seqs[get_chromosome_name(track.chrom)] + track.end, 32) ;
    return kmers ;
}

void extract_inner_kmers(vector<Track> tracks, int threads) {
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < tracks.size(); i++) {
        if (tracks[i].svtype == SVTYPE_DEL) {
            extract_deletion_kmers(tracks[i]) ;
        }
        if (tracks[i].svtype == SVTYPE_INS) {
            extract_insertion_kmers(tracks[i]) ;
        }
    }
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int main(int argc, char** argv) {
    cout << "Nebula, ultra-efficient, mapping-free genotyper." << endl ;
    if (strcmp(argv[1], "preprocess") == 0) {
        cout << "Preprocessing.." << endl ;
        KMER_TYPE = LOCUS_TYPE_INNER ;
        auto bed_tracks = load_tracks_from_file(argv[2]) ;
        string reference_path(argv[3]) ;
        int threads = std::stoi(string(argv[4]), nullptr, 10) ;
        load_chromosomes(reference_path) ;
        extract_inner_kmers(bed_tracks, threads) ;
    }
}

//int main(int argc, char** argv) {
//    string path(argv[2]) ;
//    string reference_path(argv[1]) ;
//    int threads = std::stoi(string(argv[4]), nullptr, 10) ;
//    string mode(argv[3]) ;
//    KMER_TYPE = mode == "junction" ? LOCUS_TYPE_JUNCTION : LOCUS_TYPE_INNER ;
//    time_t t ;
//    time(&t) ;
//    load_kmers(path, threads) ;
//    load_chromosomes(reference_path) ;
//    scan_reference(threads) ;
//    filter_kmers() ;
//    output_kmers(path) ;
//    time_t s ;
//    time(&s) ;
//    auto d = s - t ;
//    cout << "Returning to CgcCounterJob.." << endl ;
//}
