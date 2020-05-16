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

const int HEADER_LINE = 0 ;
const int SEQUENCE_LINE = 1 ;
const int THIRD_LINE = 2 ;
const int QUALITY_LINE = 3 ;

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

uint32_t encode_half_mer(const char* c) {
    uint32_t kmer = 0 ;
    for (uint32_t i = 0 ; i < 16 ; i++) {
        uint32_t d = (uint32_t)((c[i] & MASK) >> 1) ;
        kmer += (d << (30 - (i * 2))) ;
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

void process_read(char* seq) {
    int l = strlen(seq) ;
    l -- ; //newline skip
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
        /*
        cout << i << " " << l << endl ;
        string d_l = decode_kmer(left) ;
        string d_r = decode_kmer(right) ;
        string d_k = decode_kmer(k) ;
        cout << "===================" << endl ;
        cout << seq ;
        cout << std::setw(i) << std::right << std::setfill(' ') << std::setw(i) << (i < 32 ? "" : d_r) ;
        cout << d_k << (i + 32 + 31 < l ? d_l : "") << endl ;
        */
        if (kmer != counts->end()) {
            //cout << i << " " << l << endl ;
            //string d_l = decode_kmer(left) ;
            //string d_r = decode_kmer(right) ;
            //string d_k = decode_kmer(k) ;
            //cout << "===================" << endl ;
            //cout << seq ;
            //cout << std::setw(i) << std::right << std::setfill(' ') << std::setw(i) << (i < 32 ? "" : d_r) ;
            //cout << d_k << (i + 32 + 31 < l ? d_l : "") << endl ;
            int* total = totals->at(kmer->first) ;
            *total += 1 ;
            if (masks->find(k) == masks->end()) {
                int* count = counts->at(kmer->first) ;
                *count += 1 ;
            } else {
                std::vector<uint64_t>* m = masks->at(k) ;
                for (std::vector<uint64_t>::iterator it = m->begin(); it != m->end(); it++) {
                    //cout << decode_kmer(*it) << endl ;
                    if (i + 32 + 31 < l) {
                        if (is_subsequence(*it, left)) {
                            int* count = counts->at(kmer->first) ;
                            *count += 1 ;
                            //cout << "found" << endl ;
                            break ;
                        }
                    }
                    if (i >= 32) {
                        if (is_subsequence(*it, right)) {
                            int* count = counts->at(kmer->first) ;
                            *count += 1 ;
                            //cout << "found" << endl ;
                            break ;
                        }
                    }
                }
            }
        }
    }
}

const int SEARCHING = 0 ;
const int SKIPPING = 1 ;
const int READING = 2 ;
const unsigned long BLOCK_SIZE = 4096 * 20 ;

void output_counts(string path, int index) {
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

int get_number_of_alignments(samFile* bam_file, bam_hdr_t* bam_header, bam1_t* alignment) {
    cout << "Estimating number of reads" << endl ;
    int n = 0 ;
    time_t t ;
    time(&t) ;
    while (sam_read1(bam_file, bam_header, alignment) > 0) {
        n++ ;
    }
    time_t s ;
    time(&s) ;
    cout << "Found " << n << " reads in " << s - t << " seconds." << endl ;
    return n ;
}

// ================================ FASTQ Files ================================= \\

unordered_map<uint64_t, int> process_reads_bam(vector<bam1_t*> alignments) {
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
            /*
            cout << i << " " << l << endl ;
            string d_l = decode_kmer(left) ;
            string d_r = decode_kmer(right) ;
            string d_k = decode_kmer(k) ;
            cout << "===================" << endl ;
            cout << seq ;
            cout << std::setw(i) << std::right << std::setfill(' ') << std::setw(i) << (i < 32 ? "" : d_r) ;
            cout << d_k << (i + 32 + 31 < l ? d_l : "") << endl ;
            */
            if (kmer != counts->end()) {
                //cout << i << " " << l << endl ;
                //string d_l = decode_kmer(left) ;
                //string d_r = decode_kmer(right) ;
                //string d_k = decode_kmer(k) ;
                //cout << "===================" << endl ;
                //cout << seq ;
                //cout << std::setw(i) << std::right << std::setfill(' ') << std::setw(i) << (i < 32 ? "" : d_r) ;
                //cout << d_k << (i + 32 + 31 < l ? d_l : "") << endl ;
                int* total = totals->at(kmer->first) ;
                *total += 1 ;
                if (masks->find(k) == masks->end()) {
                    if (_counts.find(k) != _counts.end()) {
                        _counts[k] = 0 ;
                    }
                    _counts[k] += 1 ;
                    //int* count = counts->at(kmer->first) ;
                    //*count += 1 ;
                } else {
                    std::vector<uint64_t>* m = masks->at(k) ;
                    for (std::vector<uint64_t>::iterator it = m->begin(); it != m->end(); it++) {
                        //cout << decode_kmer(*it) << endl ;
                        if (i + 32 + 31 < l) {
                            if (is_subsequence(*it, left)) {
                                if (_counts.find(k) != _counts.end()) {
                                    _counts[k] = 0 ;
                                }
                                _counts[k] += 1 ;
                                //int* count = counts->at(kmer->first) ;
                                //*count += 1 ;
                                //cout << "found" << endl ;
                                break ;
                            }
                        }
                        if (i >= 32) {
                            if (is_subsequence(*it, right)) {
                                if (_counts.find(k) != _counts.end()) {
                                    _counts[k] = 0 ;
                                }
                                _counts[k] += 1 ;
                                //int* count = counts->at(kmer->first) ;
                                //*count += 1 ;
                                //cout << "found" << endl ;
                                break ;
                            }
                        }
                    }
                }
            }
        }
    }
    return _counts ;
}

vector<unordered_map<uint64_t, int>> p_process_bam(const vector<vector<bam1_t*>> &entries, int threads) {
    vector<unordered_map<uint64_t, int>> batches (threads) ;
    #pragma omp parallel for num_threads(threads)
    for(int i = 0; i <= threads; i++) {
        if (i == 0) {
            
        } else {
            unordered_map<uint64_t, int> _counts = process_reads_bam(entries[i]) ;
            batches[i] = _counts ;
        }
    }
    return batches ;
}

int process_bam(string bam, string path, int index, int threads) {
    samFile *bam_file = hts_open(bam.c_str(), "r") ;
    bam_hdr_t *bam_header = sam_hdr_read(bam_file) ; //read header
    bam1_t *alignment = bam_init1() ; //initialize an alignment
    // Timing
    time_t t ;
    time(&t) ;
    int batch_size = 1000000 ;
    int b = 0 ;
    int n = 0 ;
    uint64_t u = 0 ;
    vector<vector<bam1_t*>> entries (threads) ;
    while (sam_read1(bam_file, bam_header, alignment) > 0){
        entries[n % threads].push_back(alignment) ;
        n += 1 ;
        u += 1 ;
        if (n == batch_size) {
            n = 0 ;
            vector<unordered_map<uint64_t, int>> output = p_process_bam(entries, threads);
            for (const auto &batch : output) {
                for (const auto kmer : batch) {
                    *counts->at(kmer.first) += kmer.second ;
                }
            }
            b += 1 ;
            for (int i = 0 ; i < threads ; i++) {
                entries[i].clear() ;
            }
            time_t s ;
            time(&s) ;
            cerr << "Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << b * batch_size << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << endl ;
        }
    }
    if(!entries.empty()) {
        vector<unordered_map<uint64_t, int>> output = p_process_bam(entries, threads);
        for (const auto &batch : output) {
            for (const auto kmer : batch) {
                *counts->at(kmer.first) += kmer.second ;
            }
        }
    }
    bam_destroy1(alignment) ;
    sam_close(bam_file) ;
    output_counts(path, index) ;
    return u ;
}

// ================================ FASTQ Files ================================= \\

unordered_map<uint64_t, int> process_reads_fastq(vector<string> fastq_entries) {
    int l ;
    const char* seq ;
    unordered_map<uint64_t, int> _counts ;
    for (const auto fastq_entry : fastq_entries) {
        //TODO check these bounds
        seq = fastq_entry.c_str() ; 
        l = strlen(seq) ;
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
                //int* total = totals->at(kmer->first) ;
                //*total += 1 ;
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
    return _counts ;
}

gzFile fastq_file ;
kseq_t* fastq_iterator ;
vector<vector<vector<string>>> fastq_entries ;

bool load_next_batch(int threads, int batch_size, int p) {
    //cout << "Loading into block " << p << endl ;
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

//vector<unordered_map<uint64_t, int>> p_process_fastq(const vector<vector<string>> &entries, int threads) {
int p_process_fastq(int threads) {
    // initialize fastq_entries
    for(int i = 0; i < 2; i++) {
        fastq_entries.push_back(vector<vector<string>>(threads)) ;
    }
    int p = 0 ;
    int batch_size = 2000000 ;
    vector<unordered_map<uint64_t, int>> batches(threads) ;
    load_next_batch(threads, batch_size, p) ;
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;
    while (true) {
        // Load the next batch while processing the current one
        for (int i = 0 ; i < threads ; i++) {
            u += fastq_entries[p][i].size() ;
        }
        #pragma omp parallel for num_threads(threads + 1)
        for(int i = 0; i <= threads; i++) {
            if (i == 0) {
                load_next_batch(threads, batch_size, (p + 1) % 2) ;
            } else {
                unordered_map<uint64_t, int> _counts = process_reads_fastq(fastq_entries[p][i - 1]) ;
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

int process_fastq(string fastq, string path, int index, int threads) {
    p_process_fastq(threads) ;
    //
    //time_t t ;
    //time(&t) ;
    //int batch_size = 2000000 ;
    //int b = 0 ;
    //int l = 0 ;
    //int n = 0 ;
    //uint64_t u = 0 ;
    //while ((l = kseq_read(seq)) >= 0) {
    //    entries[n % threads].push_back(seq->seq.s) ;
    //    n += 1 ;
    //    u += 1 ;
    //    if (n == batch_size) {
    //        n = 0 ;
    //        vector<unordered_map<uint64_t, int>> output = p_process_fastq(entries, threads) ;
    //        for (const auto &batch : output) {
    //            for (const auto kmer : batch) {
    //                *counts->at(kmer.first) += kmer.second ;
    //            }
    //        }
    //        b += 1 ;
    //        for (int i = 0 ; i < threads ; i++) {
    //            entries[i].clear() ;
    //        }
    //        time_t s ;
    //        time(&s) ;
    //        cerr << "Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << b * batch_size << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << endl ;
    //    }
    //}
    //if(!entries.empty()) {
    //    vector<unordered_map<uint64_t, int>> output = p_process_fastq(entries, threads);
    //    for (const auto &batch : output) {
    //        for (const auto kmer : batch) {
    //            *counts->at(kmer.first) += kmer.second ;
    //        }
    //    }
    //}
    output_counts(path, index) ;
    return 0 ;
}


int load_kmers(int index, string path) {
    // load json file
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

int main(int argc, char** argv) {
    string path(argv[2]) ;
    string fastq(argv[3]) ;
    int index = std::stoi(string(argv[1]), nullptr, 10) ;
    int threads = std::stoi(string(argv[4]), nullptr, 10) ;
    time_t t ;
    time(&t) ;
    load_kmers(index, path) ;
    int u = 0;
    if (fastq.compare(fastq.size() - 4, 4, ".bam") == 0) {
        cout << "Input is BAM." << endl ;
        u = process_bam(fastq, path, index, threads) ;
    } else {
        cout << "Input is FASTQ." << endl ;
        fastq_file = gzopen(fastq.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
        u = process_fastq(fastq, path, index, threads) ;
    }
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    cout << "Processed " << u << " reads in " << (s - t) << " seconds. average " << u / d << " reads per second." << endl ;
}
