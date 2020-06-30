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

#include "kmer_utils.hpp"

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

vector<vector<unordered_map<uint64_t, int>>> partial_counts ;
vector<vector<unordered_map<uint64_t, int>>> partial_totals ;

void process_read(const char* seq, int l, unordered_map<uint64_t, int>& _counts, unordered_map<uint64_t, int>& _totals) {
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
            if (_totals.find(k) == _totals.end()) {
                _totals[k] = 0 ;
                _counts[k] = 0 ;
            }
            _totals[k] += 1 ;
            if (masks.find(k) == masks.end()) {
                _counts[k] += 1 ;
            } else {
                std::vector<uint64_t>* m = masks[k] ;
                for (std::vector<uint64_t>::iterator it = m->begin(); it != m->end(); it++) {
                    if (i + 32 + 31 < l) {
                        if (is_subsequence(*it, left)) {
                            _counts[k] += 1 ;
                            break ;
                        }
                    }
                    if (i >= 32) {
                        if (is_subsequence(*it, right)) {
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
// ================================ FASTQ Files ================================ \\
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

void process_batch_fastq(vector<string> fastq_entries, int thread, int p) {
    int l ;
    const char* seq ;
    for (const auto fastq_entry : fastq_entries) {
        seq = fastq_entry.c_str() ; 
        l = strlen(seq) ;
        process_read(seq, l, partial_counts[p][thread], partial_totals[p][thread]) ;
    }
}

// ============================================================================= \\
// ================================= BAM Files ================================= \\
// ============================================================================= \\

samFile *bam_file ;
bam_hdr_t *bam_header ;
bam1_t *bam_iterator ;

vector<vector<vector<bam1_t*>>> bam_entries ;

bool load_batch_bam(int threads, int batch_size, int p) {
    int i = 0 ;
    int n = 0 ;
    while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) > 0) {
        n += 1 ;
        if (n % threads == threads - 1) {
            i += 1 ;
        }
        if (n == batch_size) {
            return true ;
        }
    }
    if(bam_entries[p].empty()) {
        return false ;
    }
    return true ;
}

void process_batch_bam(vector<bam1_t*> alignments, int thread, int p) {
    char* seq = (char*) malloc(200) ;
    uint32_t len = 0 ;
    bam1_t* alignment ; 
    for (int i = 0; i < alignments.size(); i++) {
        alignment = alignments[i] ;
        uint32_t l = alignment->core.l_qseq ; //length of the read
        if (l > len) {
            if (len > 0) {
                free(seq) ;
            }
            len = l ;
            seq = (char*) malloc(l + 1) ;
        }
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++){
            seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        seq[l] = '\0' ; // null terminate
        process_read(seq, l, partial_counts[p][thread], partial_totals[p][thread]) ;
    }
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

#define BAM 0
#define FASTQ 1

int process_reads(int threads, int mode, string input_file) {
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(threads)) ;
        fastq_entries.push_back(vector<vector<string>>(threads)) ;
        partial_counts.push_back(vector<unordered_map<uint64_t, int>>(threads)) ;
        partial_totals.push_back(vector<unordered_map<uint64_t, int>>(threads)) ;
    }
    int p = 0 ;
    int batch_size = 1000000 ;
    if (mode == BAM) {
        cout << "Allocating buffers.." << endl ;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < threads; j++) {
                for (int k = 0; k <= batch_size / threads; k++) {
                    bam_entries[i][j].push_back(bam_init1()) ;
                }
            }
        }
        load_batch_bam(threads, batch_size, p) ;
    } else {
        load_batch_fastq(threads, batch_size, p) ;
    }
    cout << "Loaded initial batch.." << endl ;
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;
    while (true) {
        uint64_t v = u ;
        for (int i = 0 ; i < threads ; i++) {
            if (mode == BAM) {
                u += bam_entries[p][i].size() ;
            } else {
                u += fastq_entries[p][i].size() ;
            }
        }
        if (v == u) {
            break ;
        }
        #pragma omp parallel for num_threads(threads + 2)
        for(int i = 0; i < threads + 2; i++) {
            if (i == 0) {
                // load next batch
                if (mode == BAM) {
                    load_batch_bam(threads, batch_size, (p + 1) % 2) ;
                } else {
                    load_batch_fastq(threads, batch_size, (p + 1) % 2) ;
                }
                //cout << "Loaded." << endl ;
            } else if (i == 1) {
                // merge previous patch
                if (b >= 1) {
                    for (auto &batch : partial_counts[(p + 1) % 2]) { 
                        for (auto &kmer : batch) {
                            counts[kmer.first] += kmer.second ;
                        }
                        batch.clear() ;
                    }
                    for (auto &batch : partial_totals[(p + 1) % 2]) { 
                        for (auto &kmer : batch) {
                            totals[kmer.first] += kmer.second ;
                        }
                        batch.clear() ;
                    }
                }
            } else {
                // process current one
                if (mode == BAM) {
                    process_batch_bam(bam_entries[p][i - 2], i - 2, p) ;
                } else {
                    process_batch_fastq(fastq_entries[p][i - 2], i - 2, p) ;
                }
                //cout << "Done." << endl ;
            }
        }
        p += 1 ;
        p %= 2 ;
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
    cout << "Dumping kmer counts..." << endl ;
    string p = path + "/counts.json" ;
    cout << path << endl ;
    std::ofstream o(p);
    for (const auto kmer : counts) {
        o << decode_kmer(kmer.first) << ":" << kmer.second << ":" << totals[kmer.first] << "\n" ;
    }
}

int load_kmers(string path) {
    cout << "Loading kmers from " << path << ".." << endl ;
    ifstream json_file(path + "/pre_inner_kmers.json") ;
    nlohmann::json kmers_json ;
    json_file >> kmers_json ;
    for (nlohmann::json::iterator it = kmers_json.begin(); it != kmers_json.end(); ++it) {
        auto kmer = it.value() ;
        uint64_t k = encode_kmer(std::string(it.key()).c_str()) ;
        uint64_t rc_k = encoded_reverse_complement(k) ; 
        counts[k] = 0 ;
        counts[rc_k] = 0 ;
        totals[k] = 0 ;
        totals[rc_k] = 0 ;
        if (kmer["loci"].size()) {
            // keep this as a pointer for memory optimization, avoids keeping two copies
            std::vector<uint64_t> *m = new std::vector<uint64_t> ;
            for (nlohmann::json::iterator locus = kmer["loci"].begin(); locus != kmer["loci"].end(); ++locus) {
                for (nlohmann::json::iterator mask = locus.value()["masks"].begin(); mask != locus.value()["masks"].end(); ++mask) {
                    if (encode_kmer(mask.key()) != 0) {
                        m->push_back(encode_kmer(mask.key().c_str())) ;
                        m->push_back(encode_kmer(reverse_complement(mask.key()).c_str())) ;
                    } else {
                        //cout << "Ignoring mask.." << endl ;
                    }
                }
            }
            if (m->size()) {
                masks[k] = m ;
                masks[rc_k] = m ;
            } 
        }
    }
    cout << "Loaded " << counts.size() / 2 << " kmers." << endl ;
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
    //load_kmers(path) ;
    int u = 0;
    if (input_file.compare(input_file.size() - 4, 4, ".bam") == 0) {
        cout << "Input is BAM." << endl ;
        bam_file = hts_open(input_file.c_str(), "r") ;
        bam_header = sam_hdr_read(bam_file) ; //read header
        bam_iterator = bam_init1() ; //initialize an alignment
        u = process_reads(threads, BAM, input_file) ;
    } else {
        cout << "Input is FASTQ." << endl ;
        fastq_file = gzopen(input_file.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
        u = process_reads(threads, FASTQ, input_file) ;
    }
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    cout << "\33[2K" << "Processed " << u << " reads in " << (s - t) << " seconds. Averaging " << u / d << " reads per second." << endl ;
    output_counts(path) ;
    cout << "Returning to CgcCounterJob.." << endl ;
}
