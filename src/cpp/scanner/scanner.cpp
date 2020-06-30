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
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <pthread.h>
#include <unordered_map>

#include "kmer_utils.hpp"

#include "json.hpp"

using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

#define MASK 6

std::unordered_map<uint64_t, std::vector<uint64_t>> counts ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

#define LOCUS_TYPE_JUNCTION uint8_t(0)
#define LOCUS_TYPE_INNER    uint8_t(1)
#define LOCUS_TYPE_REF      uint8_t(2)

#define SVTYPE_MISC         uint8_t(0) 
#define SVTYPE_DEL          uint8_t(1)
#define SVTYPE_INS          uint8_t(2)

uint8_t get_chromosome_index(string chrom) {
    string ch = chrom.substr(3, chrom.length() - 3) ;
    if (ch == "X") {
        return 22 ;
    }
    if (ch == "Y") {
        return 23 ;
    }
    return std::stoi(ch) ;
}

string get_chromosome_name(uint8_t index) {
    if (index < 22) {
        return "chr" + std::to_string(index) ;
    } else if (index == 22) {
        return "chrX" ;
    } else if (index == 23) {
        return "chrY" ;
    } else {
        return "chrUn" ;
    }
}

struct Track {
    uint8_t svtype ;
    string chrom ;
    uint32_t begin ;
    uint32_t end ;
} ;

struct Locus {
    uint8_t chrom ;
    uint32_t position ;
    uint8_t type ; 
    uint64_t left ;
    uint64_t right ;
} ;

struct Kmer {
    uint32_t count ;
    uint32_t reference ;
    std::unordered_map<string, uint8_t> tracks ;
    std::vector<Locus> loci ;
    std::vector<Locus> junction_loci ;
    bool inverse ;
    Kmer(): inverse(false), count(0), reference(0) {}
} ;

void parse_locus_name(string name, Locus &locus) {
    stringstream ss(name) ;
    vector<string> tokens ;
    string token ;
    while (getline(ss, token, '_')) {
        tokens.push_back(token) ;
    }
    locus.type = LOCUS_TYPE_JUNCTION ;
    locus.chrom = get_chromosome_index(tokens[1].substr(tokens[1].find("@") + 1)) ;
    locus.position = std::stoi(tokens[2]) ;
}

Track parse_track_name(string name) {
    stringstream ss(name) ;
    vector<string> tokens ;
    string token ;
    while (getline(ss, token, '_')) {
        tokens.push_back(token) ;
    }
    int i = tokens[0].find('@') ;
    Track track ;
    track.chrom = tokens[0].substr(i, tokens[0].length() - (i + 1)) ;
    track.begin = stoi(tokens[1]) ;
    track.end = stoi(tokens[2]) ;
    return track ;
}

string get_locus_name(Locus locus) {
    if (locus.type == LOCUS_TYPE_JUNCTION) {
        return "junction_" + get_chromosome_name(locus.chrom) + "_" + std::to_string(locus.position) ;
    }
    return get_chromosome_name(locus.chrom) + "_" + std::to_string(locus.position) ;
}

std::unordered_map<string, Locus> cast_loci_map(vector<Locus>& _loci) {
    unordered_map<string, Locus> loci ;
    for (auto locus = _loci.begin(); locus != _loci.end(); locus++) {
        loci[get_locus_name(*locus)] = *locus ;
        //loci.insert(std::make_pair<string, Locus>(get_locus_name(*locus), locus)) ;
    }
    return loci ;
}

void to_json(nlohmann::json& j, const Kmer& k) {
    unordered_map<string, Locus> loci ;
    for (auto locus = k.loci.begin(); locus != k.loci.end(); locus++) {
        loci[get_locus_name(*locus)] = *locus ;
    }
    unordered_map<string, Locus> junction_loci ;
    for (auto locus = k.junction_loci.begin(); locus != k.junction_loci.end(); locus++) {
        junction_loci[get_locus_name(*locus)] = *locus ;
    }
    j = nlohmann::json{{"count", k.count}, {"reference", k.reference}, {"loci", loci}, {"junction_loci", junction_loci}, {"tracks", k.tracks}} ;
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
        //cout << locus.key() << endl ;
        Locus l = locus.value().get<Locus>() ;
        parse_locus_name(locus.key(), l) ;
        k.loci.push_back(l) ;
        //k.loci.push_back(locus.value().get<Locus>()) ;
    }
}

std::unordered_map<uint64_t, Kmer> kmers ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void process_read(string chrom, const char* seq, unordered_map<uint64_t, Kmer>& _kmers) {
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    cout << "Scanning " << chrom << ".." << endl ;
    uint64_t l = strlen(seq) ;
    for (uint64_t i = 0 ; i < l - 64 ; i++) {
        if (i == 0) {
            k = encode_kmer(seq) ;
            left = k ;
            right = encode_substring(seq, 32, 32) ;
        } else {
            k = k << 2 ;
            k += (seq[i + 31] & MASK) >> 1 ;
        }
        if (i + 32 + 31 < l) {
            right = right << 2 ;
            right += (seq[i + 32 + 31] & MASK) >> 1 ;
        }
        if (i > 32) {
            left = left << 2 ;
            left += (seq[i - 1] & MASK) >> 1 ;
        }
        if (kmers.find(k) != kmers.end()) {
            if (_kmers.find(k) == _kmers.end()) {
                _kmers[k] = Kmer() ;
                _kmers[k].count = 0 ;
            }
            string name = chrom + "@" + to_string(i) ;
            Locus locus({get_chromosome_index(chrom), uint32_t(i), LOCUS_TYPE_REF, left, right}) ;
            //_kmers[k].loci[name] = locus ;
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
                Locus locus({get_chromosome_index(chrom), uint32_t(i), LOCUS_TYPE_REF, left, right}) ;
                //_kmers[rc].loci[name] = locus ;
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
    //std::ofstream o("counts_cpp.txt") ;
    //for (auto it = kmers.begin(); it != kmers.end(); it++) {
    //    o << decode_kmer(it->first) << ":" << it->second.count << endl ;
    //}
}

vector<uint64_t> find_interest_masks(Kmer& k) {
    std::vector<uint64_t> interset_masks ;
    for (auto locus = k.loci.begin(); locus != k.loci.end(); locus++) {
        if (locus->type == LOCUS_TYPE_JUNCTION) {
            if (locus->left != 0) {
                //interset_masks.push_back(locus->second.left) ;
                interset_masks.push_back(locus->left) ;
            }
            if (locus->right != 0) {
                //interset_masks.push_back(locus->second.right) ;
                interset_masks.push_back(locus->right) ;
            }
        }
    }
    return interset_masks ;
}

bool is_kmer_returning(Kmer& kmer) {
    for(auto track = kmer.tracks.begin(); track != kmer.tracks.end(); track++) {
        Track t = parse_track_name(track->first) ;
        //cout << t.chrom << "_" << t.begin << "_" << t.end << endl ;
        for(auto locus = kmer.loci.begin(); locus != kmer.loci.end(); locus++) {
            if (locus->type != LOCUS_TYPE_JUNCTION) {
                if (abs(int(locus->position) - int(t.end)) < 250) {
                    return true ;
                }
                if (abs(int(locus->position) - int(t.begin)) < 250) {
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
    // actual filtering
    // Parallelize this
    //cout << "Filtering " << kmers.size() << " kmers based on masks.." << endl ;
    kmer = kmers.begin() ;
    while (kmer != kmers.end()) {
        Kmer& k = kmer->second ;
        std::vector<uint64_t> interset_masks = find_interest_masks(k) ;
        //std::unordered_map<string, Locus> loci(k.loci) ;
        std::vector<Locus> loci(k.loci) ;
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
                locus = k.loci.erase(locus) ;
                continue ;
            }
            locus++ ;
        }
        int l_2 = k.loci.size() ;
        // loci with less than two masks exist
        //auto a = find_if(k.loci.begin(), k.loci.end(), [](std::pair<string, Locus> l) {
        auto a = find_if(k.loci.begin(), k.loci.end(), [](Locus l) {
            return l.left == 0 || l.right == 0 ;
        }) ;
        //non-junction loci exists
        //auto b = find_if(k.loci.begin(), k.loci.end(), [](std::pair<string, Locus> l) {
        auto b = find_if(k.loci.begin(), k.loci.end(), [](Locus l) {
            return l.type != LOCUS_TYPE_JUNCTION ;
        }) ;
        // loci exist wtih less than two masks
        if (a != k.loci.end()) {
            if (l_1 != l_2) { // some non-junction loci were filtered
                if (b == k.loci.end()) { // all non-junction loci were filtered
                    k.loci = loci ;
                    auto locus = k.loci.begin() ; // count non-junction loci only
                    while (locus != k.loci.end()) {
                        if (locus->type == LOCUS_TYPE_JUNCTION) {
                            k.junction_loci.push_back(*locus) ;
                            locus = k.loci.erase(locus) ;
                        } else {
                            locus++ ;
                        }
                    }
                    // TODO: set inverse to True
                } else { // smome non-junction loci remain
                    kmer = kmers.erase(kmer) ;
                    continue ;
                }
            } else {
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
    string p = path + "/FilterJunctionKmersJob/kmers.json" ;
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
// ================================ FASTA Files ================================ \\
// ============================================================================= \\

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
                             //reuse this line
                             break ;
                         }
                         for (int i = 0; i < line.length(); i++) {
                            line[i] = toupper(line[i]) ;
                         }
                         //cout << line << endl ;
                         //std::this_thread::sleep_for(std::chrono::seconds(1)) ;
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
                    // don't read line again
                    continue ;
                }
            }
        } 
        if (!std::getline(fasta_file, line)) {
            break;
        }
    }
    free(buffer) ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int parse_svtype(string svtype) {
    if (svtype == "DEL") {
        return SVTYPE_DEL ;
    }
    if (svtype == "INS") {
        return SVTYPE_INS ;
    }
    return SVTYPE_MISC ;
}

int load_kmers(string path, int threads) {
    cout << "Loading kmers from " << path << ".." << endl ;
    ifstream json_file(path + "/ExtractJunctionKmersJob/batch_merge.json") ;
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
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < l; i++) {
        ifstream track_file(filenames[i]) ;
        nlohmann::json track_json ;
        track_file >> track_json ;
        int t = omp_get_thread_num() ;
        for (nlohmann::json::iterator kmer = track_json.begin(); kmer != track_json.end(); kmer++) {
            uint64_t k = encode_kmer(canonicalize(kmer.key()).c_str()) ;
            _kmers[t][k] = kmer.value().get<Kmer>() ;
            _kmers[t][k].count = 0 ;
            _kmers[t][k].tracks[tracknames[i]] = 1 ;
        }
        track_file.close() ;
        n += 1 ;
        if (t == 0) {
            if (n % 100 == 0) {
                cout << "Loaded " << threads << "x" << n << " tracks approximately.." << '\r' << flush ;
            }
        }
    }
    cout << endl ;
    json_file.close() ;
    for (int i = 0; i < threads; i++) {
        for (auto kmer = _kmers[i].begin(); kmer != _kmers[i].end(); kmer++) {
            if (kmers.find(kmer->first) == kmers.end()) {
                kmers[kmer->first] = kmer->second ;
            }
        }
    }
    cout << "Loaded " << kmers.size() << " kmers." << endl ;
    return 0 ;
}
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int main(int argc, char** argv) {
    string path(argv[2]) ;
    string reference_path(argv[1]) ;
    int threads = std::stoi(string(argv[3]), nullptr, 10) ;
    time_t t ;
    time(&t) ;
    load_kmers(path, threads) ;
    load_chromosomes(reference_path) ;
    scan_reference(threads) ;
    filter_kmers() ;
    output_kmers(path) ;
    int u = 0;
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    cout << "Returning to CgcCounterJob.." << endl ;
}
