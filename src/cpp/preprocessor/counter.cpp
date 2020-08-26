#include <omp.h>
#include <ctime>
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

#include "counter.hpp"

#include "json.hpp"

using namespace std ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void KmerCounter::process_read(const char* seq, int l, unordered_map<uint64_t, int>& _counts, unordered_map<uint64_t, int>& _totals) {
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    //cout << seq << endl ;
    for (int i = 0 ; i <= l - 32 ; i++) {
        if (i == 0) {
            k = encode_kmer(seq) ;
            right = encode_kmer(seq + 32) ;
            left = k ;
        } else {
            k = k << 2 ;
            k += (seq[i + 31] & MASK) >> 1 ;
        }
        if (i > 32) {
            left = left << 2 ;
            left += (seq[i - 1] & MASK) >> 1 ;
        }
        if (i + 32 + 31 < l) {
            right = right << 2 ;
            right += (seq[i + 32 + 31] & MASK) >> 1 ;
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
                    if (i >= 32) {
                        if (is_subsequence(*it, left)) {
                            _counts[k] += 1 ;
                            break ;
                        }
                    }
                    if (i + 32 + 31 < l) {
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

bool KmerCounter::load_batch_fastq(int p) {
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

void KmerCounter::process_batch_fastq(vector<string> fastq_entries, int thread, int p) {
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

bool KmerCounter::load_batch_bam(int p) {
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
    return n != 0 ? true : false ;
}

void KmerCounter::process_batch_bam(vector<bam1_t*> alignments, int thread, int p) {
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

void KmerCounter::process_reads(string input_file) {
    auto c = Configuration::getInstance() ;
    if (input_file.compare(input_file.size() - 4, 4, ".bam") == 0) {
        cout << "Input is BAM." << endl ;
        mode = COUNTER_MODE_BAM ;
        bam_file = hts_open(input_file.c_str(), "r") ;
        bam_header = sam_hdr_read(bam_file) ; //read header
    } else {
        cout << "Input is FASTQ." << endl ;
        mode = COUNTER_MODE_FASTQ ;
        fastq_file = gzopen(input_file.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
        batch_size = 1000000 ;
    }
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(threads)) ;
        fastq_entries.push_back(vector<vector<string>>(threads)) ;
        partial_counts.push_back(vector<unordered_map<uint64_t, int>>(threads)) ;
        partial_totals.push_back(vector<unordered_map<uint64_t, int>>(threads)) ;
        batch_size = 1000000 ;
    }
    int p = 0 ;
    if (mode == COUNTER_MODE_BAM) {
        cout << "Allocating buffers.." << endl ;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < threads; j++) {
                for (int k = 0; k <= batch_size / threads; k++) {
                    bam_entries[i][j].push_back(bam_init1()) ;
                }
            }
        }
        load_batch_bam(p) ;
    } else {
        load_batch_fastq(p) ;
    }
    cout << "Loaded initial batch.." << endl ;
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;
    bool should_continue = true ;
    while (true) {
        uint64_t v = u ;
        for (int i = 0; i < threads; i++) {
            if (mode == COUNTER_MODE_BAM) {
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
                if (mode == COUNTER_MODE_BAM) {
                    should_continue = load_batch_bam((p + 1) % 2) ;
                } else {
                    load_batch_fastq((p + 1) % 2) ;
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
                if (mode == COUNTER_MODE_BAM) {
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
        if (!should_continue) {
            break ;
        }
    }
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    cout << "\33[2K" << "Processed " << u << " reads in " << (s - t) << " seconds. Averaging " << u / d << " reads per second." << endl ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void KmerCounter::run() {
    auto c = Configuration::getInstance() ;
    load_kmers() ;
    if (c->bam.size() != 0) {
        for (auto it = c->bam.begin(); it != c->bam.end(); it++) {
            process_reads(*it) ;
        }
    }
    if (c->fastq.size() != 0) {
        for (auto it = c->fastq.begin(); it != c->fastq.end(); it++) {
            process_reads(*it) ;
        }
    }
    output_counts() ;
}

void KmerCounter::load_kmers() {
    vector<unordered_map<uint64_t, Kmer>*> kmers{depth_kmers, gc_kmers, genotyping_kmers} ;
    for (auto _kmers = kmers.begin(); _kmers != kmers.end(); _kmers++) {
        //cout << counts.size() << endl ;
        for (auto it = (*_kmers)->begin(); it != (*_kmers)->end(); it++) {
            uint64_t k = it->first ;
            uint64_t rc_k = encoded_reverse_complement(k) ; 
            counts[k] = 0 ;
            counts[rc_k] = 0 ;
            totals[k] = 0 ;
            totals[rc_k] = 0 ;
            auto kmer = it->second ;
            if (kmer.loci.size()) {
                std::vector<uint64_t> *m = new std::vector<uint64_t> ;
                for (auto locus = kmer.loci.begin(); locus != kmer.loci.end(); ++locus) {
                    if (locus->left != 0) {
                        m->push_back(locus->left) ;
                        m->push_back(encoded_reverse_complement(locus->left)) ;
                    }
                    if (locus->right != 0) {
                        m->push_back(locus->right) ;
                        m->push_back(encoded_reverse_complement(locus->right)) ;
                    }
                }
                if (m->size()) {
                    masks[k] = m ;
                    masks[rc_k] = m ;
                } 
            }
        }
    }
    cout << "Loaded " << counts.size() / 2 << " kmers with " << masks.size() << " masks." << endl ;
}

void KmerCounter::output_counts() {
    auto c = Configuration::getInstance() ;
    string p = c->workdir + "/counts.txt" ;
    cout << "Dumping kmer counts to " << p << ".." << endl ;
    std::ofstream o(p);
    for (const auto kmer : counts) {
        o << decode_kmer(kmer.first) << ":" << kmer.second << ":" << totals[kmer.first] << "\n" ;
    }
    nlohmann::json payload ;
    for (auto kmer: counts) {
        uint64_t canon = encode_kmer(canonicalize(decode_kmer(kmer.first))) ;
        if (gc_kmers->find(canon) != gc_kmers->end()) {
            gc_kmers->find(canon)->second.count += kmer.second ;
        }
        if (depth_kmers->find(canon) != depth_kmers->end()) {
            depth_kmers->find(canon)->second.count += kmer.second ;
        }
        if (genotyping_kmers->find(canon) != genotyping_kmers->end()) {
            genotyping_kmers->find(canon)->second.count += kmer.second ;
        }
    }
}

void KmerCounter::load_counts() {
    auto c = Configuration::getInstance() ;
    cout << "Loading kmer counts.." << endl ;
    int num_batches = 1000 ;
    #pragma omp parallel for num_threads(c->threads)
    for (int i = 0; i < num_batches; i++) {
        string path = c->workdir + "/counts_batch_" + std::to_string(i) ;
        std::ifstream o(path) ;
        string line ;
        while (std::getline(o, line)) {
            stringstream ss(line) ;
            vector<string> tokens ;
            string token ;
            while (getline(ss, token, ':')) {
                tokens.push_back(token) ;
            }
            int count = std::stoi(tokens[1]) ;
            int total = std::stoi(tokens[2]) ;
            uint64_t canon = encode_kmer(canonicalize(tokens[0])) ;
            if (gc_kmers->find(canon) != gc_kmers->end()) {
                gc_kmers->find(canon)->second.count += count ;
            }
            if (depth_kmers->find(canon) != depth_kmers->end()) {
                depth_kmers->find(canon)->second.count += count ;
            }
            if (genotyping_kmers->find(canon) != genotyping_kmers->end()) {
                genotyping_kmers->find(canon)->second.count += count ;
                genotyping_kmers->find(canon)->second.total += total ;
            }
        }
    }
}
