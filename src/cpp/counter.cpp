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

void handler(int sig) {
    void *array[10];
    size_t size;
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

#define MASK 6

std::unordered_map<uint64_t, std::unordered_map<std::string, string>*> *half_mers = new std::unordered_map<uint64_t, std::unordered_map<std::string, string>*> ;
std::unordered_map<uint64_t, bool> *other_mers = new std::unordered_map<uint64_t, bool> ;
std::unordered_map<uint64_t, int*> *counts = new std::unordered_map<uint64_t, int*> ;
std::unordered_map<uint64_t, int*> *gaps = new std::unordered_map<uint64_t, int*> ;
std::unordered_map<uint64_t, std::vector<uint64_t>*> *masks = new std::unordered_map<uint64_t, std::vector<uint64_t>*> ;

//string reverse_complement(string s) {
//    std::replace(s.begin(), s.end(), 'A', 'Z') ;
//    std::replace(s.begin(), s.end(), 'T', 'A') ;
//    std::replace(s.begin(), s.end(), 'Z', 'T') ;
//    std::replace(s.begin(), s.end(), 'C', 'Z') ;
//    std::replace(s.begin(), s.end(), 'G', 'C') ;
//    std::replace(s.begin(), s.end(), 'Z', 'G') ;
//    std::reverse(s.begin(), s.end());
//    return s ;
//}

//uint8_t encode_base(char base) {
//    if (base == 'A') {
//        return (uint8_t)0 ;
//    }
//    if (base == 'T') {
//        return (uint8_t)1 ;
//    }
//    if (base == 'C') {
//        return (uint8_t)2 ;
//    }
//    return (uint8_t)3 ;
//}
//
//// Will encode 4 bases in a byte, 32 bases in 8 bytes = 64bit uint64_t
//void encode_read(char* read) {
//    int l = strlen(read) ;
//    int m ;
//    if (l % 4 == 0){
//        m = l / 4 ;
//    } else {
//        m = l / 4 + 1 ;
//    }
//    uint8_t* a = (uint8_t*) malloc(sizeof(uint8_t) * m) ;
//    int r = 0 ;
//    int j = 0 ;
//    a[0] = 0 ;
//    for (int i = 0 ; i < l - 1 ; i++) {
//        uint8_t b = encode_base(read[i]) ;
//        a[j] = a[j] & (b << (3 - j) * 2) ;
//        r += 1 ;
//        if (r == 4) {
//            j += 1 ;
//            r = 0 ;
//            a[j] = 0 ;
//        }
//    }
//}

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
        else {
            cout << "BAD KMER" << endl ;
            cout << "#" << s << "#" << endl ;
            int j ;
            cin >> j ;
        }
        j += 1 ;
    }
    return std::string(rc) ;
}

string canonicalize(string s) {
    string rc = reverse_complement(s) ;
    return s.compare(rc) <= 0 ? s : rc ;
}

bool is_subsequence(std::string* x, std::string* y) {
    string::iterator x_it = x->begin();
    for (string::iterator y_it = y->begin(); y_it != y->end() && x_it != x->end(); ++y_it) {
        if (*x_it == *y_it) {
            ++x_it ;
        }
    }
    bool res = (x_it == x->end()) ;
    return res ;
}

uint64_t encode_kmer(const char* c) {
    uint64_t kmer = 0 ;
    for (int i = 0 ; i < 32 ; i++) {
        //uint64_t r = 0 ? c[i] == 'A' : 1 ? c[i] == 'C' : 2 ? c[i] == 'G' : 3 ;
        //kmer += (r << (62 - (i << 1))) ;
        kmer += ((uint64_t)((c[i] & MASK) >> 1) << (62 - (i << 1))) ;
    }
    return kmer ;
}

uint64_t encode_substring(const char* c, int base, int l) {
    uint64_t mask = 0 ;
    for (int i = 0 ; i < l ; i++) {
        //uint64_t r = 0 ? c[base + i] == 'A' : 1 ? c[base + i] == 'C' : 2 ? c[base + i] == 'G' : 3 ;
        //mask += (r << (62 - (i << 1))) ;
        mask += ((uint64_t)((c[base + i] & MASK) >> 1) << (62 - (i << 1))) ;
    }
    return mask ;
}

string decode_kmer(uint64_t kmer) {
    char* d = (char*) malloc(sizeof(char) * 32) ;
    uint64_t mask = 3 ;
    for (int i = 0; i < 32 ; i++) {
        uint64_t t = kmer & mask ;
        if (t == 1) {
            d[31 - i] = 'C' ;
        }
        else if (t == 0) {
            d[31 - i] = 'A' ;
        }
        else if (t == 2) {
            d[31 - i] = 'T' ;
        }
        else if (t == 3) {
            d[31 - i] = 'G' ;
        }
        else {
            cout << "error" << endl ;
        }
        kmer = kmer >> 2 ;
    }
    cout << d << endl ;
    return std::string(d) ;
}

void process_read(char* read, int index) {
    int l = strlen(read) ;
    char c ;
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    for (int i = 0 ; i <= l - 32 ; i++) {
        if (i == 0) {
            k = encode_kmer(read) ;
            left = encode_substring(read, 32, 32) ;
            right = k ;
        } else {
            k = k << 2 ;
            k += (read[i + 31] & MASK) >> 1 ;
        }
        if (i + 32 + 31 < l) {
            left = left << 2 ;
            left += (read[i + 32 + 31] & MASK) >> 1 ;
        }
        if (i > 32) {
            right = right << 2 ;
            right += (read[i - 1] & MASK) >> 1 ;
        }
        auto kmer = counts->find(k) ;
        /*
        cout << i << endl ;
        string d_l = decode_kmer(left) ;
        string d_r = decode_kmer(right) ;
        string d_k = decode_kmer(k) ;
        ios init(NULL);
        init.copyfmt(cout);
        cout << "===================" << endl ;
        cout << read ;
        cout << std::setw(i) << std::right << std::setfill(' ') << std::setw(i) << (i < 32 ? "" : d_r) ;
        cout.copyfmt(init) ;
        cout << d_k << (i + 32 + 31 < l ? d_l : "") << endl ;
        */
        if (kmer != counts->end()) {
            std::vector<uint64_t>* m = masks->at(k) ;
            for (std::vector<uint64_t>::iterator it = m->begin(); it != m->end(); it++) {
                if (i + 32 + 31 < l) {
                    if (__builtin_popcountll((left ^ *it)) <= 3) {
                        int* count = counts->at(kmer->first) ;
                        *count += 1 ;
                        break ;
                    }
                }
                if (i >= 32) {
                    if (__builtin_popcountll((right ^ *it)) <= 3) {
                        int* count = counts->at(kmer->first) ;
                        *count += 1 ;
                        break ;
                    }
                }
            }
        }
    }
}

//void process_gapped_read(char* read, int num, bool select) {
//    int l = strlen(read) ;
//    std::string seq(read) ;
//    std::unordered_map<std::string, std::vector<int>> index ;
//    for (int i = 0 ; i <= l - 1 - 15 ; i++) {
//        std::string k = seq.substr(i, 15) ;
//        if (half_mers->find(k) == half_mers->end() and other_mers->find(k) == other_mers->end()) {
//            continue ;
//        }
//        index[k].push_back(i) ;
//    }
//    for (std::unordered_map<string, std::vector<int>>::iterator it = index.begin(); it != index.end(); it++) {
//        auto half = half_mers->find(it->first) ;
//        if (half == half_mers->end()) {
//            continue ;
//        }
//        for(std::unordered_map<std::string, std::string>::iterator j = half->second->begin(); j != half->second->end(); j++){
//            std::string kmer = j->second ;
//            auto other = index.find(j->first) ;
//            if (other != index.end()) {
//                for (std::vector<int>::iterator a = it->second.begin(); a != it->second.end(); a++) {
//                    for (std::vector<int>::iterator b = other->second.begin(); b != other->second.end(); b++) {
//                        int d = *b - (*a + 15) ;
//                        if (select) {
//                            if (d >= 0 && d <= 10) {
//                                int* count = counts->at(kmer) ;
//                                count[d] = count[d] + 1 ;
//                            }
//                            continue ;
//                        }
//                        if (d >= 0 && d == *gaps->at(kmer)) {
//                            int* count = counts->at(kmer) ;
//                            *count += 1 ;
//                        }
//                    }
//                }
//            }
//        }
//    }
//}

int process_fastq(string fastq, string path, int index, int threads, int type) {
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
                if (s - t != 0) {
                    cout << std::left << setw(2) << index << " progress: " << setw(14) << std::fixed << p ;
                    cout << " took: " << setw(7) << std::fixed << s - t << " ETA: " << setw(14) << e ;
                    cout << " current: " << setw(12) << ftell(fastq_file) << " limit: " << (index + 1) * chunk_size ;
                    cout << " reads per second: " << u / (s - t) << " Done" << endl ;
                }
            }
            if (type == 0) {
                process_read(line, u) ;
            } else {
                //process_gapped_read(line, u, type == 2) ;
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
        if (r == -1) {
            break ;
        }
    }
    //json payload ;
    //string p = path + "/c_batch_" + std::to_string(index) + ".json" ;
    //std::ofstream o(p);
    //for (std::unordered_map<string, int*>::iterator i = counts->begin(); i != counts->end(); i++) {
    //    int* count = i->second ;
    //    if (type == 2) {
    //        int c = -1 ;
    //        int g ;
    //        for (int j = 0 ; j < 10 ; j++) {
    //            if (count[j] >= c) {
    //                g = j ;
    //                c = count[j] ;
    //            }
    //        }
    //        o << std::string(i->first) << ":" << c << ":" << g << endl ;
    //    } else {
    //        o << std::string(i->first) << ":" << *count << endl ;
    //    }
    //}
    return u ;
}

int transform(int index, string path, string fastq, int threads, int type) {
    cout << "index " << index << " counting " << fastq << endl ;
    ifstream json_file(path + "/pre_kmers.json") ;
    json kmers_json ;
    json_file >> kmers_json ;
    for (json::iterator it = kmers_json.begin(); it != kmers_json.end(); ++it) {
        auto kmer = it.value() ;
        uint64_t k = encode_kmer(std::string(it.key()).c_str()) ;
        uint64_t rc_k = encode_kmer(reverse_complement(std::string(it.key())).c_str()) ;
        int* count = new int ;
        *count = 0 ;
        counts->emplace(std::make_pair(k, count)) ;
        counts->emplace(std::make_pair(rc_k, count)) ;
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
    cout << counts->size() / 2 << " kmers" << endl ;
    return 0 ;
}

//int transform_gapped(int index, string path, string fastq, int threads, int type) {
//    cout << "index " << index << " counting " << fastq << endl ;
//    cout << path << "/half_mers.json" << endl ;
//    // load json file
//    cout << "loading kmers" << endl ;
//    ifstream kmers_json_file(path + "/pre_kmers.json") ;
//    json kmers_json ;
//    kmers_json_file >> kmers_json ;
//    //
//    for (json::iterator i = kmers_json.begin(); i != kmers_json.end(); ++i) {
//        auto kmer = i.value() ;
//        if (type == 2) {
//            int* count = new int[11] ;
//            counts->insert(std::make_pair(std::string(i.key()), count)) ;
//        } else {
//            int* count = new int ;
//            *count = 0 ;
//            counts->insert(std::make_pair(std::string(i.key()), count)) ;
//            int* gap = new int ;
//            *gap = std::stoi(kmer["gap"].get<std::string>(), nullptr, 10) ;
//            gaps->insert(std::make_pair(std::string(i.key()), gap)) ;
//        }
//    }
//    // 
//    cout << "loading half mers" << endl ;
//    ifstream half_mer_json_file(path + "/half_mers.json") ;
//    json half_mers_json ;
//    half_mer_json_file >> half_mers_json ;
//    //
//    for (json::iterator i = half_mers_json.begin(); i != half_mers_json.end(); ++i) {
//        auto half_mer = i.value() ;
//        std::unordered_map<string, string>* map = new std::unordered_map<string, string> ;
//        for (json::iterator j = half_mer.begin(); j != half_mer.end(); ++j) { 
//            map->insert(std::make_pair(std::string(j.key()), std::string(j.value().get<std::string>()))) ;
//            other_mers->insert(std::make_pair(std::string(j.key()), true)) ;
//        }
//        half_mers->insert(std::make_pair(std::string(i.key()), map)) ;
//    }
//    cout << counts->size() << " kmers" << endl ;
//    return 0 ;
//} 

int main(int argc, char** argv) {
    cout << ('C' & MASK) << endl ;
    cout << ('A' & MASK) << endl ;
    cout << ('T' & MASK) << endl ;
    cout << ('G' & MASK) << endl ;
    string path(argv[2]) ;
    string fastq(argv[3]) ;
    int index = std::stoi(string(argv[1]), nullptr, 10) ;
    int threads = std::stoi(string(argv[4]), nullptr, 10) ;
    int type = std::stoi(string(argv[5]), nullptr, 10) ;
    cout << fastq << endl ;
    cout << "Threads: " << threads << endl ;
    time_t t ;
    time(&t) ;
    int n ;
    if (type == 0) {
        transform(index, path, fastq, threads, type) ;
    } else {
        //transform_gapped(index, path, fastq, threads, type) ;
    }
    n = process_fastq(fastq, path, index, threads, type) ;
    time_t s;
    time(&s) ;
    auto d = s - t ;
    cout << "processed " << n << " reads in " << (s - t) << " seconds. average" << n / d << " reads per second." << endl ;
}
