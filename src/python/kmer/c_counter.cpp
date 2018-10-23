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

// reduction
std::unordered_map<std::string, std::unordered_map<std::string, void*>*> *kmers = new std::unordered_map<std::string, std::unordered_map<std::string, void*>*> ;
// gapped
std::unordered_map<std::string, std::unordered_map<std::string, string>*> *half_mers = new std::unordered_map<std::string, std::unordered_map<std::string, string>*> ;
std::unordered_map<std::string, bool> *other_mers = new std::unordered_map<std::string, bool> ;
std::unordered_map<std::string, int*> *counts = new std::unordered_map<std::string, int*> ;
std::unordered_map<std::string, int*> *gaps = new std::unordered_map<std::string, int*> ;
std::unordered_map<std::string, std::vector<string>*> *masks = new std::unordered_map<std::string, std::vector<string>*> ;

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
        if (s[i] == 'T') {
            rc[j] = 'A' ;
        }
        if (s[i] == 'C') {
            rc[j] = 'G' ;
        }
        if (s[i] == 'G') {
            rc[j] = 'C' ;
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

//bool is_canonical_subsequence(string* x, string* y) {
//    return is_subsequence(reverse_complement(x), y) || is_subsequence(x, y) ;
//}

void process_read(char* read, int index) {
    int l = strlen(read) ;
    std::string seq(read) ;
    int slack = (100 - 31) / 2 ;
    // for each kmer in the read
    for (int i = 0 ; i <= l - 1 - 31 ; i++) {
        std::string k = seq.substr(i, 31) ;
        auto kmer = kmers->find(k) ;
        // check if kmer exists
        if (kmer != kmers->end()) {
            std::vector<string>* m = static_cast<std::vector<string>*>(kmer->second->at("masks")) ;
            std::string right = seq.substr(i + 31, slack) ;
            std::string left = i <= slack ? seq.substr(0, i) : seq.substr(i - slack, slack) ;
            for (std::vector<string>::iterator it = m->begin(); it != m->end(); it++) {
                if (is_subsequence(&(*it), &left) || is_subsequence(&(*it), &right)) {
                    int* count = counts->at(kmer->first) ;
                    *count += 1 ;
                    //cout << *it << endl ;
                    //cout << k << endl ;
                    //cout << left << endl ;
                    //cout << right << endl ;
                    //cout << index << endl ;
                    //cout << "waiting for input..." << endl ;
                    //int j ;
                    //cin >> j ;
                    break ;
                }
            }
        }
    }
}

void process_gapped_read(char* read, int num) {
    if (num != 686) {
        return ;
    }
    cout << num << endl ;
    // get rid of the EOL
    int l = strlen(read) ;
    std::string seq(read) ;
    // index kmers
    cout << seq << endl ;
    std::unordered_map<std::string, std::vector<int>> index ;
    for (int i = 0 ; i <= l - 1 - 15 ; i++) {
        cout << 1 << endl ;
        std::string k = seq.substr(i, 15) ;
        cout << k << endl ;
        if (half_mers->find(k) == half_mers->end() and other_mers->find(k) == other_mers->end()) {
            cout << "not found" << endl ;
            continue ;
        }
        if (index.find(k) == index.end()) {
            cout << "adding to index" << endl ;
            index.emplace(std::make_pair(k, std::vector<int>(1, i))) ;
        } else {
            cout << 6 << endl ;
            index[k].push_back(i) ;
        }
    }
    cout << index.size() << endl ;
    // for every half-mer in index
    for (std::unordered_map<string, std::vector<int>>::iterator it = index.begin(); it != index.end(); it++) {
        cout << 5 << endl ;
        auto half = half_mers->find(it->first) ;
        cout << 2 << endl ;
        if (half == half_mers->end()) {
            continue ;
        }
        cout << 3 << endl ;
        //for all second pairs of this half
        for(std::unordered_map<std::string, std::string>::iterator j = half->second->begin(); j != half->second->end(); j++){
            std::string kmer = j->second ;
            auto other = index.find(j->first) ;
            if (other != index.end()) {
                for (std::vector<int>::iterator a = it->second.begin(); a != it->second.end(); a++) {
                    for (std::vector<int>::iterator b = other->second.begin(); b != other->second.end(); b++) {
                        cout << 4 << endl ;
                        int d = *b - (*a + 15) ;
                        if (d == *gaps->at(kmer)) {
                            //cout << "count" << endl ;
                            int* count = counts->at(kmer) ;
                            *count += 1 ;
                            //cout << num << endl ;
                            //cout << it->first << endl ;
                            //cout << j->first << endl ;
                            //cout << "waiting for input..." << endl ;
                            //int y ;
                            //cin >> y ;
                            //cout << it->first << endl ;
                            //cout << j->first << endl ;
                            //cout << kmer << endl ;
                            //cout << num << endl ;
                            //cout << "waiting for input..." << endl ;
                            //int j ;
                            //cin >> j ;
                            //std::string right = seq.substr(*b + 15, 25) ;
                            //std::string left = *a <= 25 ? seq.substr(0, *a) : seq.substr(*a - 25, 25) ;
                            //auto m = masks->at(kmer) ;
                            //for (std::vector<string>::iterator it = m->begin(); it != m->end(); it++) {
                            //    if (is_subsequence(&(*it), &left) || is_subsequence(&(*it), &right)) {
                            //        break ;
                            //    }
                            //}
                        }
                    }
                }
            }
        }
    }
}


void process_fastq(string fastq, string path, int index, int threads, int type) {
    std::ifstream in(fastq, std::ifstream::ate | std::ifstream::binary) ;
    unsigned long size = (unsigned long) in.tellg() ;
    cout << "Size: " << size << endl ;
    unsigned long chunk_size = size / 100 ;
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
                time_t s;
                time(&s) ;
                double p = c / double(chunk_size) ;
                double e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600) ;
                cout.precision(10) ;
                cout << std::left << setw(2) << index << " progress: " << setw(14) << std::fixed << p ;
                cout << " took: " << setw(7) << std::fixed << s - t << " ETA: " << setw(14) << e ;
                cout << " current: " << setw(12) << ftell(fastq_file) << " limit: " << (index + 1) * chunk_size << endl ;
            }
            if (type == 0) {
                process_read(line, u) ;
            } else {
                process_gapped_read(line, u) ;
            }
        }
        else if (state == THIRD_LINE) {
            state = QUALITY_LINE ;
        }
        else if (state == QUALITY_LINE) {
            state = HEADER_LINE ;
        }
        line = strcpy(line, ahead) ;
        r = getline(&ahead, &len_ahead, fastq_file) ;
        if (r == -1) {
            break ;
        }
    }
    json payload ;
    string p = path + "/c_batch_" + std::to_string(index) + ".json" ;
    std::ofstream o(p);
    for (std::unordered_map<string, int*>::iterator i = counts->begin(); i != counts->end(); i++) {
        int* count = i->second ;
        o << std::string(i->first) << ":" << *count << endl ;
    }
}

int transform(int index, string path, string fastq, int threads, int type) {
    cout << "index " << index << " counting " << fastq << endl ;
    ifstream json_file(path + "/pre_kmers.json") ;
    json kmers_json ;
    json_file >> kmers_json ;
    for (json::iterator it = kmers_json.begin(); it != kmers_json.end(); ++it) {
        auto kmer = it.value() ;
        std::unordered_map<string, void*>* map = new std::unordered_map<string, void*> ;
        int* count = new int ;
        *count = 0 ;
        int* doubt = new int ;
        *doubt = 0 ;
        map->emplace(std::make_pair("count", count)) ;
        map->emplace(std::make_pair("doubt", doubt)) ;
        counts->emplace(std::make_pair(std::string(it.key()), count)) ;
        counts->emplace(std::make_pair(reverse_complement(std::string(it.key())), count)) ;
        std::vector<std::string> *m = new std::vector<std::string> ;
        for (json::iterator locus = kmer["loci"].begin(); locus != kmer["loci"].end(); ++locus) {
            for (json::iterator mask = locus.value()["masks"].begin(); mask != locus.value()["masks"].end(); ++mask) {
                m->push_back(mask.key()) ;
                m->push_back(reverse_complement(mask.key())) ;
            }
        }
        map->emplace(std::make_pair("masks", m)) ;
        kmers->emplace(std::make_pair(std::string(it.key()), map)) ;
        kmers->emplace(std::make_pair(reverse_complement(std::string(it.key())), map)) ;
    }
    cout << kmers->size() << " kmers" << endl ;
    process_fastq(fastq, path, index, threads, type) ;
}

void transform_gapped(int index, string path, string fastq, int threads, int type) {
    cout << "index " << index << " counting " << fastq << endl ;
    cout << path << "/half_mers.json" << endl ;
    // load kmer masks
    ifstream kmers_json_file(path + "/pre_kmers.json") ;
    json kmers_json ;
    kmers_json_file >> kmers_json ;
    for (json::iterator i = kmers_json.begin(); i != kmers_json.end(); ++i) {
        auto kmer = i.value() ;
        std::vector<std::string> *m = new std::vector<std::string> ;
        for (json::iterator mask = kmer["masks"].begin(); mask != kmer["masks"].end(); ++mask) {
            m->push_back(mask.key()) ;
            m->push_back(reverse_complement(mask.key())) ;
        }
        masks->emplace(std::make_pair(std::string(i.key()), m)) ;
        int* count = new int ;
        *count = 0 ;
        counts->emplace(std::make_pair(std::string(i.key()), count)) ;
        int* gap = new int ;
        *gap = std::stoi(kmer["gap"].get<std::string>(), nullptr, 10) ;
        gaps->emplace(std::make_pair(std::string(i.key()), gap)) ;
    }
    cout << masks->size() << " kmers" << endl ;
    // load half-mer pairs
    ifstream half_mer_json_file(path + "/half_mers.json") ;
    json half_mers_json ;
    half_mer_json_file >> half_mers_json ;
    for (json::iterator i = half_mers_json.begin(); i != half_mers_json.end(); ++i) {
        auto half_mer = i.value() ;
        std::unordered_map<string, string>* map = new std::unordered_map<string, string> ;
        for (json::iterator j = half_mer.begin(); j != half_mer.end(); ++j) { 
            map->emplace(std::make_pair(std::string(j.key()), std::string(j.value().get<std::string>()))) ;
            other_mers->emplace(std::make_pair(std::string(j.key()), true)) ;
        }
        half_mers->emplace(std::make_pair(std::string(i.key()), map)) ;
    }
    cout << counts->size() << " kmers" << endl ;
    //
    process_fastq(fastq, path, index, threads, type) ;
} 

int main(int argc, char** argv) {
    //signal(SIGSEGV, handler);
    string path(argv[2]) ;
    string fastq(argv[3]) ;
    int index = std::stoi(string(argv[1]), nullptr, 10) ;
    int threads = std::stoi(string(argv[4]), nullptr, 10) ;
    int type = std::stoi(string(argv[5]), nullptr, 10) ;
    cout << fastq << endl ;
    cout << "Threads: " << threads << endl ;
    if (type == 0) {
        transform(index, path, fastq, threads, type) ;
    } else {
        transform_gapped(index, path, fastq, threads, type) ;
    }
}
