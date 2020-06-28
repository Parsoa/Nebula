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
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <hts_internal.h>
#include "htslib/hts_endian.h"

#include "htslib/hfile.h"

#include "kseq.h"
#include "json.hpp"

KSEQ_INIT(gzFile, gzread)
using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

#define MASK 6

std::unordered_map<uint64_t, int> counts ;
std::unordered_map<uint64_t, int> totals ;
std::unordered_map<uint64_t, std::vector<uint64_t>*> masks ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\


// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

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
    int m = 0 ;
    int n = 0 ;
    int j = 46 ;
    uint64_t mask_y = MASK_Y ;
    uint64_t mask_x = MASK_X ;
    x = x >> 8 ;
    for (int i = 62; i >= 0 && j >= 0; i -= 2) {
        if ((x & mask_x) >> j == (y & mask_y) >> i) {
            j -= 2 ;
            mask_x = mask_x >> 2 ;
        } else {
            if (n < 2) {
                if ((x & (mask_x >> 2)) >> (j - 2) == (y & (mask_y >> 2)) >> (i - 2)) {
                    j -= 2 ;
                    mask_x = mask_x >> 2 ;
                }
            }
        }
        mask_y = mask_y >> 2 ; 
    }
    bool res = j == -2 ;
    return res ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void process_read(const char* seq, int l, unordered_map<uint64_t, int>& _counts) {
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    //cout << seq << endl ;
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
        if (counts.find(k) != counts.end()) {
            totals[k] += 1 ;
            if (masks.find(k) == masks.end()) {
                if (_counts.find(k) == _counts.end()) {
                    _counts[k] = 0 ;
                }
                _counts[k] += 1 ;
            } else {
                std::vector<uint64_t>* m = masks[k] ;
                for (std::vector<uint64_t>::iterator it = m->begin(); it != m->end(); it++) {
                    if (i + 32 + 31 < l) {
                        if (is_subsequence(*it, left)) {
                            if (_counts.find(k) == _counts.end()) {
                                _counts[k] = 0 ;
                            }
                            _counts[k] += 1 ;
                            break ;
                        }
                    }
                    if (i >= 32) {
                        if (is_subsequence(*it, right)) {
                            if (_counts.find(k) == _counts.end()) {
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
// ================================ FASTA Files ================================ \\
// ============================================================================= \\

unordered_map<string, char*> chromosomes ;

int get_reference_size(ifstream f) {
    fasta_file.seekg(0, ios_base::end) ;
    int l = fasta_file.tellg() ;
    fasta_file.seekg(0, ios_base::begin) ;
    return l ;
}

void load_chromosomes(string path) {
    // assume human genome length
    ifstream fasta_file ;
    fasta_file.open(path, ios::binary) ;
    int l = get_reference_size(fasta_file) ;
    char* reference = (char*) malloc(sizeof(char) * l) ;
    //
    char* buffer = (char*) malloc(sizeof(char) * 4096) ;
    // read all of file
    std::string line ;
    while (std::getline(fasta_file, line)) {
        if (line.substr(0, 4) == ">chr") {
            int n = line.length() ;
            if (n == 
        }
    }
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
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    output_counts(path) ;
    cout << "Returning to CgcCounterJob.." << endl ;
}
