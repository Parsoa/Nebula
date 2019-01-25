#include <locale>
#include <ctime>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <bitset>

#include "json.hpp"

using namespace std ;
using json = nlohmann::json ;

int JOB = 0 ;
bool DEBUG = false ;
uint64_t MASK_Y = (uint64_t) 3 << 62 ;
uint64_t MASK_X = (uint64_t) 3 << 46 ;

#define MASK 6

int SELECT_GAPPED_KMERS = 2 ;
int COUNT_GAPPED_KMERS = 1 ;
int COUNT_INNER_KMERS = 0 ;
int COUNT_KMERS = 3 ;

std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint64_t>*> *half_mers = new std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint64_t>*> ;
std::unordered_map<uint32_t, bool> *other_mers = new std::unordered_map<uint32_t, bool> ;
std::unordered_map<uint64_t, int*> *counts = new std::unordered_map<uint64_t, int*> ;
std::unordered_map<uint64_t, int*> *totals = new std::unordered_map<uint64_t, int*> ;
std::unordered_map<uint64_t, int*> *gaps = new std::unordered_map<uint64_t, int*> ;
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

bool is_subsequence(uint64_t x, uint64_t y) {
    x = x >> 8 ;
    //cout << "^ " << decode_kmer(y) << endl ;
    //cout << "# " << decode_kmer(x) << endl ;
    uint64_t mask_y = MASK_Y ;
    uint64_t mask_x = MASK_X ;
    int j = 46 ;
    int m = 0 ;
    int n = 0 ;
    int u = 0 ;
    for (int i = 62; i >= 0 && j >= 0; i -= 2) {
        if ((x & mask_x) >> j == (y & mask_y) >> i) {
            //cout << "(^" << decode_base((y & mask_y) >> i) << "," << decode_base((x & mask_x) >> j) << ")" << endl ;
            j -= 2 ;
            mask_x = mask_x >> 2 ;
            //cout << decode_kmer(mask_x) << " " << decode_kmer((x & mask_x) >> j) << "," << j << endl ;
        } else {
            if (n < 2) {
                if ((x & (mask_x >> 2)) >> (j - 2) == (y & (mask_y >> 2)) >> (i - 2)) {
                    j -= 2 ;
                    mask_x = mask_x >> 2 ;
                }
            }
            //cout << "(#" << decode_base((y & mask_y) >> i) << "," << decode_base((x & mask_x) >> j) << ")" << endl ;
            //cout << decode_kmer(mask_x) << " " << decode_kmer((x & mask_x) >> j) << "," << j << endl ;
        }
        mask_y = mask_y >> 2 ; 
    }
    bool res = j == -2 ;
    /*if (res) {
        cout << "subsequence" << endl ;
        cin >> j ;
    } else {
        cout << "not" << endl ;
    }*/
    return res ;
}

void process_read(char* seq, int read, int line) {
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
            if (JOB == COUNT_INNER_KMERS) {
                //cout << line << " " << read << endl ;
                //cout << i << " " << l << endl ;
                //string d_l = decode_kmer(left) ;
                //string d_r = decode_kmer(right) ;
                //string d_k = decode_kmer(k) ;
                //cout << "===================" << endl ;
                //cout << seq ;
                //cout << std::setw(i) << std::right << std::setfill(' ') << std::setw(i) << (i < 32 ? "" : d_r) ;
                //cout << d_k << (i + 32 + 31 < l ? d_l : "") << endl ;
                std::vector<uint64_t>* m = masks->at(k) ;
                int* total = totals->at(kmer->first) ;
                *total += 1 ;
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
            } else {
                int* count = counts->at(kmer->first) ;
                *count += 1 ;
            }
        }
    }
}

void process_gapped_read(char* seq, int read, int line) {
    int l = strlen(seq) ;
    l-- ;
    char c ;
    uint32_t k = 0 ;
    std::unordered_map<uint32_t, std::vector<int>> index ;
    bool found_half = false ;
    for (int i = 0 ; i <= l - 16 ; i++) {
        if (i == 0) {
            k = encode_half_mer(seq) ;
        } else {
            k = k << 2 ;
            k += (seq[i + 15] & MASK) >> 1 ;
        }
        if (half_mers->find(k) == half_mers->end()) {
            if (found_half == true) {
                if (other_mers->find(k) != other_mers->end()) {
                    index[k].push_back(i) ;
                }
            }
        } else {
            // don't need to store other kmers as long as we haven't seen a valid half_mer yet.
            found_half = true ;
            index[k].push_back(i) ;
        }
    }
    //
    for (std::unordered_map<uint32_t, std::vector<int>>::iterator it = index.begin(); it != index.end(); it++) {
        auto half = half_mers->find(it->first) ;
        if (half == half_mers->end()) {
            continue ;
        }
        for(std::unordered_map<uint32_t, uint64_t>::iterator j = half->second->begin(); j != half->second->end(); j++){
            uint64_t kmer = j->second ;
            auto other = index.find(j->first) ;
            if (other != index.end()) {
                for (std::vector<int>::iterator a = it->second.begin(); a != it->second.end(); a++) {
                    for (std::vector<int>::iterator b = other->second.begin(); b != other->second.end(); b++) {
                        int d = *b - (*a + 16) ;
                        if (JOB == SELECT_GAPPED_KMERS) {
                            if (d >= 0 && d <= 10) {
                                int* count = counts->at(kmer) ;
                                count[d] = count[d] + 1 ;
                            }
                        } else {
                            if (d >= 0 && d == *gaps->at(kmer)) {
                                int* count = counts->at(kmer) ;
                                *count += 1 ;
                            }
                        }
                    }
                }
            }
        }
    }
}

void output_counts(string path, int index) {
    json payload ;
    cout << "dumping kmer counts..." << endl ;
    string p = path + "/c_batch_" + std::to_string(index) + ".json" ;
    cout << path << endl ;
    std::ofstream o(p);
    for (std::unordered_map<uint64_t, int*>::iterator i = counts->begin(); i != counts->end(); i++) {
        int* count = i->second ;
        if (JOB == SELECT_GAPPED_KMERS) {
            o << decode_kmer(i->first) ;
            for (int i = 0 ; i <= 10; i++) {
                o << ":" << count[i] ;
            }
            o << endl ;
        } else if (JOB == COUNT_GAPPED_KMERS) {
            o << decode_kmer(i->first) << ":" << *count << endl ; 
        } else {
            auto total = totals->find(i->first) ;
            o << decode_kmer(i->first) << ":" << *count << ":" << *total->second << endl ;
        }
    }
    cout << "done" << endl ;
}

int process_fastq(string fastq, string path, int index, int threads) {
    std::ifstream in(fastq, std::ifstream::ate | std::ifstream::binary) ;
    unsigned long size = (unsigned long) in.tellg() ;
    cout << "Size: " << size << endl ;
    cout << "Threads:" << threads << endl ;
    unsigned long chunk_size = size / threads ;
    unsigned long offset = index * chunk_size ;
    cout << "Chunk size: " << chunk_size << endl ;
    cout << "Offset: " << offset << endl ;
    FILE* fastq_file = fopen(fastq.c_str(), "r") ;
    fseek(fastq_file, offset, SEEK_SET) ;
    const int HEADER_LINE = 0 ;
    const int SEQUENCE_LINE = 1 ;
    const int THIRD_LINE = 2 ;
    const int QUALITY_LINE = 3 ;
    std::locale loc ;
    size_t len_line = 0 ;
    size_t len_ahead = 0 ;
    char* line = NULL ;
    char* ahead = NULL ;
    int r = 0 ;
    int state = HEADER_LINE ;
    getline(&line, &len_line, fastq_file) ; 
    getline(&ahead, &len_ahead, fastq_file) ; 
    int n = 0 ;
    int m = 0 ;
    int u = 0 ;
    int v = 2 ;
    time_t t ;
    time(&t) ;
    while (true) {
        if (state == HEADER_LINE) {
            if (line[0] == '@' && ahead[0] != '@') {
                long int l = ftell(fastq_file) ;
                if (l >= (index + 1) * chunk_size) {
                    cout << "Index " << index << " reached segment boundary" << endl ;
                    break ;
                }
                state = SEQUENCE_LINE ;
            }
        }
        else if (state == SEQUENCE_LINE) {
            state = THIRD_LINE ;
            n += 1 ;
            u += 1 ;
            if (n == 100000) {
                n = 0 ;
                m += 1 ;
                unsigned long c = ftell(fastq_file) - index * chunk_size ;
                time_t s ;
                time(&s) ;
                double p = c / double(chunk_size) ;
                double e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600) ;
                cout.precision(10) ;
                if (s - t != 0 & DEBUG == 0) {
                    cout << std::left << setw(2) << index << " progress: " << setw(14) << std::fixed << p ;
                    cout << " took: " << setw(7) << std::fixed << s - t << " ETA: " << setw(14) << e ;
                    cout << " current: " << setw(12) << ftell(fastq_file) << " limit: " << (index + 1) * chunk_size ;
                    cout << " reads per second: " << u / (s - t) << endl ;
                }
            }
            if (JOB == COUNT_INNER_KMERS || JOB == COUNT_KMERS) {
                process_read(line, u, v) ;
            } else {
                process_gapped_read(line, u, v) ;
            }
        }
        else if (state == THIRD_LINE) {
            state = QUALITY_LINE ;
        }
        else if (state == QUALITY_LINE) {
            state = HEADER_LINE ;
        }
        char* tmp = line ;
        line = ahead ;
        ahead = tmp ;
        r = getline(&ahead, &len_ahead, fastq_file) ;
        v++ ;
        if (r == -1) {
            break ;
        }
    }
    if (DEBUG == 0) {
        output_counts(path, index) ;
    }
    return u ;
}

int transform(int index, string path) {
    // load json file
    cout << "loading kmers" << endl ;
    ifstream json_file(path + "/pre_kmers.json") ;
    json kmers_json ;
    json_file >> kmers_json ;
    for (json::iterator it = kmers_json.begin(); it != kmers_json.end(); ++it) {
        auto kmer = it.value() ;
        uint64_t k = encode_kmer(std::string(it.key()).c_str()) ;
        //uint64_t k = encode_kmer(std::string("GAGATTGCGCCATTGCACTCCAGCCTGGGCGT").c_str()) ;
        //uint64_t k = encode_kmer(std::string("CCTGTCTCTACTAAAAATACAGATAAATTAGC").c_str()) ;
        //cout << decode_kmer(k) << " " << std::string(it.key()) << endl ;
        //if (decode_kmer(k).compare(std::string(it.key())) != 0) {
        //    cout << "Error on kmer #" << std::string(it.key()) << "#" << decode_kmer(k) << "#" << endl ;
        //}
        uint64_t rc_k = encode_kmer(reverse_complement(std::string(it.key())).c_str()) ;
        //uint64_t rc_k = encode_kmer(reverse_complement(std::string("GAGATTGCGCCATTGCACTCCAGCCTGGGCGT")).c_str()) ;
        //uint64_t rc_k = encode_kmer(reverse_complement(std::string("CCTGTCTCTACTAAAAATACAGATAAATTAGC")).c_str()) ;
        //cout << decode_kmer(rc_k) << " " << reverse_complement(std::string(it.key())) << endl ;
        //if (decode_kmer(rc_k).compare(reverse_complement(std::string(it.key()))) != 0) {
        //    cout << "Error on kmer #" << reverse_complement(std::string(it.key())) << "#" << decode_kmer(rc_k) << "#" << endl ;
        //}
        int* count = new int ;
        *count = 0 ;
        counts->emplace(std::make_pair(k, count)) ;
        counts->emplace(std::make_pair(rc_k, count)) ;
        int* total = new int ;
        *total = 0 ;
        totals->emplace(std::make_pair(k, total)) ;
        totals->emplace(std::make_pair(rc_k, total)) ;
        k = encode_kmer(std::string(it.key()).c_str()) ;
        rc_k = encode_kmer(reverse_complement(std::string(it.key())).c_str()) ;
        if (JOB == COUNT_INNER_KMERS) {
            std::vector<uint64_t> *m = new std::vector<uint64_t> ;
            for (json::iterator locus = kmer["loci"].begin(); locus != kmer["loci"].end(); ++locus) {
                for (json::iterator mask = locus.value()["masks"].begin(); mask != locus.value()["masks"].end(); ++mask) {
                    m->push_back(encode_kmer(mask.key().c_str())) ;
                    m->push_back(encode_kmer(reverse_complement(mask.key()).c_str())) ;
                }
            }
            masks->emplace(std::make_pair(k, m)) ;
            masks->emplace(std::make_pair(rc_k, m)) ;
        }
    }
    cout << counts->size() / 2 << " kmers" << endl ;
    return 0 ;
}

int transform_gapped(int index, string path) {
    // load json file
    cout << "loading kmers" << endl ;
    ifstream kmers_json_file(path + "/pre_kmers.json") ;
    json kmers_json ;
    kmers_json_file >> kmers_json ;
    //
    for (json::iterator i = kmers_json.begin(); i != kmers_json.end(); ++i) {
        auto kmer = i.value() ;
        uint64_t k = encode_kmer(std::string(i.key()).c_str()) ;
        if (JOB == SELECT_GAPPED_KMERS) {
            int* count = new int[11] ;
            counts->insert(std::make_pair(k, count)) ;
            for (int j = 0 ; j <= 10 ; j++) {
                count[j] = 0;
            }
        } else {
            int* count = new int ;
            *count = 0 ;
            counts->insert(std::make_pair(k, count)) ;
            int* gap = new int ;
            //*gap = std::stoi(kmer["gap"].get<std::string>(), nullptr, 10) ;
            *gap = kmer["gap"] ;
            gaps->insert(std::make_pair(k, gap)) ;
        }
    }
    // 
    cout << "loading half mers" << endl ;
    ifstream half_mer_json_file(path + "/half_mers.json") ;
    json half_mers_json ;
    half_mer_json_file >> half_mers_json ;
    //
    for (json::iterator i = half_mers_json.begin(); i != half_mers_json.end(); ++i) {
        auto half_mer = i.value() ;
        std::unordered_map<uint32_t, uint64_t>* map = new std::unordered_map<uint32_t, uint64_t> ;
        for (json::iterator j = half_mer.begin(); j != half_mer.end(); ++j) { 
            uint32_t h = encode_half_mer(std::string(j.key()).c_str()) ;
            map->insert(std::make_pair(h, encode_kmer(std::string(j.value().get<std::string>()).c_str()))) ;
            other_mers->insert(std::make_pair(h, true)) ;
        }
        half_mers->insert(std::make_pair(encode_half_mer(std::string(i.key()).c_str()), map)) ;
    }
    cout << counts->size() << " kmers" << endl ;
    return 0 ;
} 

int main(int argc, char** argv) {
    //cout << ('C' & MASK) << endl ;
    //cout << ('A' & MASK) << endl ;
    //cout << ('T' & MASK) << endl ;
    //cout << ('G' & MASK) << endl ;
    //uint64_t b = encode_kmer(std::string("GCTCACTGGCAACTTCTGCCTCCCCGGTTCAA").c_str()) ;
    //uint64_t a = encode_kmer(std::string("AGCTCACTGCAACTTCTGCCTCCCCGGTTCAA").c_str()) ;
    //cout << is_subsequence(a, b) << endl ;
    string path(argv[2]) ;
    string fastq(argv[3]) ;
    int index = std::stoi(string(argv[1]), nullptr, 10) ;
    int threads = std::stoi(string(argv[4]), nullptr, 10) ;
    JOB = std::stoi(string(argv[5]), nullptr, 10) ;
    DEBUG = std::stoi(string(argv[6]), nullptr, 10) ;
    cout << fastq << endl ;
    cout << "Threads: " << threads << endl ;
    time_t t ;
    time(&t) ;
    cout << "index " << index << " counting " << fastq << endl ;
    if (JOB == COUNT_INNER_KMERS || JOB == COUNT_KMERS) {
        transform(index, path) ;
    } else if (JOB == COUNT_GAPPED_KMERS || JOB == SELECT_GAPPED_KMERS) {
        transform_gapped(index, path) ;
    }
    int n = process_fastq(fastq, path, index, threads) ;
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    cout << "processed " << n << " reads in " << (s - t) << " seconds. average " << n / d << " reads per second." << endl ;
}
