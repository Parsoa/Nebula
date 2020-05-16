#include <omp.h>
#include <ctime>
#include <string>
#include <thread>
#include <locale>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <pthread.h>
#include <unordered_map>

#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "kseq.h"
#include "json.hpp"

KSEQ_INIT(gzFile, gzread)
using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#  define NUM_THREADS 1
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#  define NUM_THREADS 16
#endif

#define MASK 6

std::unordered_map<uint64_t, int*> *counts = new std::unordered_map<uint64_t, int*> ;
std::unordered_map<uint64_t, int*> *totals = new std::unordered_map<uint64_t, int*> ;
std::unordered_map<uint64_t, std::vector<uint64_t>*> *masks = new std::unordered_map<uint64_t, std::vector<uint64_t>*> ;

string reverse_complement(string s) {
    char rc[s.length()] ;
    int j = 0 ;
    for (int i = s.length() - 1 ; i >= 0 ; i--) {
        if (s[i] == 'A') {
            rc[j] = 'T' ;
        }
        else if (s[i] == 'T') {
            rc[j] = 'A' ;
        }
        else if (s[i] == 'C') {
            rc[j] = 'G' ;
        }
        else if (s[i] == 'G') {
            rc[j] = 'C' ;
        }
        else if (s[i] == 'N') {
            rc[j] == 'N' ;
        }
        j += 1 ;
    }
    return std::string(rc) ;
}

string canonicalize(string s) {
    string rc = reverse_complement(s) ;
    return s.compare(rc) <= 0 ? s : rc ;
}

uint64_t encode_kmer(const char* c) {
    uint64_t kmer = 0 ;
    for (uint64_t i = 0 ; i < 32 ; i++) {
        uint64_t d = (uint64_t)((c[i] & MASK) >> 1) ;
        kmer += (d << (62 - (i * 2))) ;
    }
    return kmer ;
}

uint64_t encode_substring(const char* c, int base, int l) {
    uint64_t mask = 0 ;
    for (uint64_t i = 0 ; i < l ; i++) {
        mask += ((uint64_t)((c[base + i] & MASK) >> 1) << (62 - (i * 2))) ;
    }
    return mask ;
}

string decode_kmer(uint64_t kmer) {
    char* d = (char*) malloc(sizeof(char) * 33) ;
    d[32] = 0 ;
    uint64_t mask = 3 ;
    for (int i = 0; i < 32 ; i++) {
        uint64_t t = kmer & mask ;
        if (t == 0) {
            d[31 - i] = 'A' ;
        }
        else if (t == 1) {
            d[31 - i] = 'C' ;
        }
        else if (t == 2) {
            d[31 - i] = 'T' ;
        }
        else if (t == 3) {
            d[31 - i] = 'G' ;
        }
        kmer = kmer >> 2 ;
    }
    return std::string(d) ;
}

string decode_base(uint64_t base) {
    if (base == 0) {
        return "A" ;
    }
    else if (base == 1) {
        return "C" ;
    }
    else if (base == 2) {
        return "T" ;
    }
    else if (base == 3) {
        return "G" ;
    }
    return "N" ;
}

uint64_t MASK_Y = (uint64_t) 3 << 62 ;
uint64_t MASK_X = (uint64_t) 3 << 46 ;

bool is_subsequence(uint64_t x, uint64_t y) {
    x = x >> 8 ;
    DEBUG(cout << "^ " << decode_kmer(y) << endl ;)
    DEBUG(cout << "# " << decode_kmer(x) << endl ;)
    uint64_t mask_y = MASK_Y ;
    uint64_t mask_x = MASK_X ;
    int j = 46 ;
    int m = 0 ;
    int n = 0 ;
    for (int i = 62; i >= 0 && j >= 0; i -= 2) {
        if ((x & mask_x) >> j == (y & mask_y) >> i) {
            DEBUG(cout << "(^" << decode_base((y & mask_y) >> i) << "," << decode_base((x & mask_x) >> j) << ")" << endl ;)
            j -= 2 ;
            mask_x = mask_x >> 2 ;
            DEBUG(cout << decode_kmer(mask_x) << " " << decode_kmer((x & mask_x) >> j) << "," << j << endl ;)
        } else {
            DEBUG(cout << "(#" << decode_base((y & mask_y) >> i) << "," << decode_base((x & mask_x) >> j) << ")" << endl ;)
            DEBUG(cout << decode_kmer(mask_x) << " " << decode_kmer((x & mask_x) >> j) << "," << j << endl ;)
            if (n < 2) {
                if ((x & (mask_x >> 2)) >> (j - 2) == (y & (mask_y >> 2)) >> (i - 2)) {
                    DEBUG(cout << "(>" << decode_base((y & (mask_y >> 2)) >> (i - 2)) << "," << decode_base((x & (mask_x >> 2)) >> (j - 2)) << ")" << endl ;)
                    j -= 2 ;
                    mask_x = mask_x >> 2 ;
                }
            }
        }
        mask_y = mask_y >> 2 ; 
    }
    bool res = j == -2 ;
    #ifdef DEBUG_MODE
    if (res) {
        cout << "subsequence" << endl ;
        cin >> j ;
    } else {
        cout << "not" << endl ;
    }
    #endif
    return res ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

#define BAM 0
#define FASTQ 1

void process_read(const char* seq, int l, unordered_map<uint64_t, int>& _counts) {
    char c ;
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    for (int i = 0 ; i <= l - 32 ; i++) {
        if (i == 0) {
            k = encode_kmer(seq) ;
            left = encode_substring(seq, 32, 32) ;
            right = k ;
        } else {
            k = k << 2 ;
            k += (seq[i + 31] & MASK) >> 1 ;
        }
        if (i + 32 + 31 < l) {
            left = left << 2 ;
            left += (seq[i + 32 + 31] & MASK) >> 1 ;
        }
        if (i > 32) {
            right = right << 2 ;
            right += (seq[i - 1] & MASK) >> 1 ;
        }
        auto kmer = counts->find(k) ;
        if (kmer != counts->end()) {
            int* total = totals->at(kmer->first) ;
            *total += 1 ;
            if (masks->find(k) == masks->end()) {
                if (_counts.find(k) != _counts.end()) {
                    _counts[k] = 0 ;
                }
                _counts[k] += 1 ;
            } else {
                std::vector<uint64_t>* m = masks->at(k) ;
                for (std::vector<uint64_t>::iterator it = m->begin(); it != m->end(); it++) {
                    if (i + 32 + 31 < l) {
                        if (is_subsequence(*it, left)) {
                            if (_counts.find(k) != _counts.end()) {
                                _counts[k] = 0 ;
                            }
                            _counts[k] += 1 ;
                            break ;
                        }
                    }
                    if (i >= 32) {
                        if (is_subsequence(*it, right)) {
                            if (_counts.find(k) != _counts.end()) {
                                _counts[k] = 0 ;
                            }
                            _counts[k] += 1 ;
                            break ;
                        }
                    }
                }
            }
        }
    }
}

// ============================================================================= \\
// ================================= BAM Files ================================== \\
// ============================================================================= \\

samFile *bam_file ;
bam_hdr_t *bam_header ;
bam1_t *bam_iterator ;
vector<vector<vector<bam1_t*>>> bam_entries ;

bool load_batch_bam(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        bam_entries[p][i].clear() ;
    }
    int n = 0 ;
    while (sam_read1(bam_file, bam_header, bam_iterator) > 0){
        bam_entries[p][n % threads].push_back(bam_iterator) ;
        n += 1 ;
        if (n == batch_size) {
            return true ;
        }
    }
    if(bam_entries[p].empty()) {
        return false ;
    }
    return true ;
}

unordered_map<uint64_t, int> process_batch_bam(vector<bam1_t*> alignments) {
    char* seq ; //= (char*) malloc(200) ;
    uint32_t len = 0 ;
    unordered_map<uint64_t, int> _counts ;
    for (const auto alignment : alignments) {
        uint32_t l = alignment->core.l_qseq ; //length of the read
        if (l > len) {
            if (len > 0) {
                free(seq) ;
            }
            len = l ;
            seq = (char*) malloc(l + 1) ;
        }
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < len; i++){
            seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        seq[l] = '\0' ; // null terminate
        process_read(seq, l, _counts) ;
    }
    return _counts ;
}

// ============================================================================= \\
// ================================ FASTQ Files ================================= \\
// ============================================================================= \\

gzFile fastq_file ;
kseq_t* fastq_iterator ;
vector<vector<vector<string>>> fastq_entries ;

bool load_batch_fastq(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        fastq_entries[p][i].clear() ;
    }
    int l = 0 ;
    int n = 0 ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        fastq_entries[p][n % threads].push_back(fastq_iterator->seq.s) ;
        n += 1 ;
        if (n == batch_size) {
            return true ;
        }
    }
    if(fastq_entries[p].empty()) {
        return false ;
    }
    return true ;
}

unordered_map<uint64_t, int> process_batch_fastq(vector<string> fastq_entries) {
    int l ;
    const char* seq ;
    unordered_map<uint64_t, int> _counts ;
    for (const auto fastq_entry : fastq_entries) {
        seq = fastq_entry.c_str() ; 
        l = strlen(seq) ;
        process_read(seq, l, _counts) ;
    }
    return _counts ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int process_reads(int threads, int mode) {
    // initialize bam_entries
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(threads)) ;
        fastq_entries.push_back(vector<vector<string>>(threads)) ;
    }
    int p = 0 ;
    int batch_size = 2000000 ;
    vector<unordered_map<uint64_t, int>> batches(threads) ;
    if (mode == BAM) {
        load_batch_bam(threads, batch_size, p) ;
    } else {
        load_batch_fastq(threads, batch_size, p) ;
    }
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;
    while (true) {
        // Load the next batch while processing the current one
        for (int i = 0 ; i < threads ; i++) {
            if (mode == BAM) {
                u += bam_entries[p][i].size() ;
            } else {
                u += fastq_entries[p][i].size() ;
            }
        }
        #pragma omp parallel for num_threads(threads + 1)
        for(int i = 0; i <= threads; i++) {
            if (i == 0) {
                if (mode == BAM) {
                    load_batch_bam(threads, batch_size, (p + 1) % 2) ;
                } else {
                    load_batch_fastq(threads, batch_size, (p + 1) % 2) ;
                }
            } else {
                unordered_map<uint64_t, int> _counts ;
                if (mode == BAM) {
                    _counts = process_batch_bam(bam_entries[p][i - 1]) ;
                } else {
                    _counts = process_batch_fastq(fastq_entries[p][i - 1]) ;
                }
                batches[i - 1] = _counts ;
            }
        }
        p += 1 ;
        p %= 2 ;
        for (const auto &batch : batches) {
            for (const auto kmer : batch) {
                *counts->at(kmer.first) += kmer.second ;
            }
        }
        b += 1 ;
        time_t s ;
        time(&s) ;
        cerr << "Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << u << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << "\r" ;
    }
    return u ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void output_counts(string path) {
    nlohmann::json payload ;
    cout << "dumping kmer counts..." << endl ;
    string p = path + "/count.json" ;
    cout << path << endl ;
    std::ofstream o(p);
    for (auto i = counts->begin(); i != counts->end(); i++) {
        int* count = i->second ;
        auto total = totals->find(i->first) ;
        o << decode_kmer(i->first) << ":" << *count << ":" << *total->second << "\n" ;
    }
    cout << "done" << endl ;
}

int load_kmers(string path) {
    cout << "Loading kmers from " << path << ".." << endl ;
    ifstream json_file(path + "/pre_inner_kmers.json") ;
    nlohmann::json kmers_json ;
    json_file >> kmers_json ;
    for (nlohmann::json::iterator it = kmers_json.begin(); it != kmers_json.end(); ++it) {
        auto kmer = it.value() ;
        uint64_t k = encode_kmer(std::string(it.key()).c_str()) ;
        uint64_t rc_k = encode_kmer(reverse_complement(std::string(it.key())).c_str()) ;
        int* count = new int ;
        *count = 0 ;
        counts->emplace(std::make_pair(k, count)) ;
        counts->emplace(std::make_pair(rc_k, count)) ;
        int* total = new int ;
        *total = 0 ;
        totals->emplace(std::make_pair(k, total)) ;
        totals->emplace(std::make_pair(rc_k, total)) ;
        if (kmer["loci"].size()) {
            std::vector<uint64_t> *m = new std::vector<uint64_t> ;
            for (nlohmann::json::iterator locus = kmer["loci"].begin(); locus != kmer["loci"].end(); ++locus) {
                for (nlohmann::json::iterator mask = locus.value()["masks"].begin(); mask != locus.value()["masks"].end(); ++mask) {
                    m->push_back(encode_kmer(mask.key().c_str())) ;
                    m->push_back(encode_kmer(reverse_complement(mask.key()).c_str())) ;
                }
            }
            if (m->size()) {
                masks->emplace(std::make_pair(k, m)) ;
                masks->emplace(std::make_pair(rc_k, m)) ;
            } 
        }
    }
    cout << "Loaded " << counts->size() / 2 << " kmers." << endl ;
    return 0 ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int main(int argc, char** argv) {
    string path(argv[1]) ;
    string input_file(argv[2]) ;
    int threads = std::stoi(string(argv[3]), nullptr, 10) ;
    time_t t ;
    time(&t) ;
    load_kmers(path) ;
    int u = 0;
    if (input_file.compare(input_file.size() - 4, 4, ".bam") == 0) {
        cout << "Input is BAM." << endl ;
        bam_file = hts_open(input_file.c_str(), "r") ;
        bam_header = sam_hdr_read(bam_file) ; //read header
        bam_iterator = bam_init1() ; //initialize an alignment
        u = process_reads(threads, BAM) ;
    } else {
        cout << "Input is FASTQ." << endl ;
        fastq_file = gzopen(input_file.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
        u = process_reads(threads, FASTQ) ;
    }
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    output_counts(path) ;
    cout << "Processed " << u << " reads in " << (s - t) << " seconds. average " << u / d << " reads per second." << endl ;
}
