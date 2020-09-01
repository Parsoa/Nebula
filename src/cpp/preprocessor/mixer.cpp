#include <omp.h>
#include <mutex>
#include <math.h>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_map>

#include "mixer.hpp"
#include "config.hpp"

#include "json.hpp"

using namespace std ;

void Mixer::run() {
    auto c = Configuration::getInstance() ;
    load_tracks() ;
    int n = 0 ;
    int m = 0 ;
    for (auto track = bed_tracks.begin(); track != bed_tracks.end(); track++) {
        kmers[*track] ;
        flags[*track] = false ;
    }
    #pragma omp parallel for num_threads(c->threads)
    for (int i = 0; i < bed_tracks.size(); i++) {
        bool any_positive = false ;
        auto track = bed_tracks[i] ;
        for (auto it = tracks.begin(); it != tracks.end(); it++) {
            if (it->find(track) != it->end()) {
                if (it->find(track)->first.genotype.find("1") != string::npos) {
                    any_positive = true ;
                }
            }
        }
        if (any_positive) {
            m += 1 ;
            load_kmers(track) ;
        }
        n += 1 ;
        if (n % 100 == 0) {
            cout << "Progress " << std::left << std::setw(6) << std::fixed << std::setprecision(3) << float(n) / bed_tracks.size() << "%..\r" << flush ;
        }
    }
    cout << endl ;
    cout << n << endl ;
    cout << m << endl ;
    cout << "--------------------------------------------------------- " << endl ;
    export_kmers() ;
}

void Mixer::load_tracks() {
    auto c = Configuration::getInstance() ;
    for (auto it = c->bed.begin(); it != c->bed.end(); it++) {
        tracks.push_back(load_tracks_from_file_as_dict(*it)) ;
    }
    for(auto sample: c->samples) {
        string path = c->workdir + "/" + sample + "/genotypes.bed" ;
        calls.push_back(load_tracks_from_file_as_dict(path)) ;
    }
    int n = 0 ;
    int m = 0 ;
    int p = 0 ;
    int q = 0 ;
    for (auto track = bed_tracks.begin(); track != bed_tracks.end(); track++) {
        if (calls[0].find(*track) != calls[0].end()) {
            n++ ;
        }
        if (calls[1].find(*track) != calls[1].end()) {
            m++ ;
        }
        if (tracks[0].find(*track) != tracks[0].end()) {
            p++ ;
        }
        if (tracks[1].find(*track) != tracks[1].end()) {
            q++ ;
        }
    }
    cout << n << " in HG00514" << endl ;
    cout << m << " in HG00733" << endl ;
    cout << p << " in HG00514" << endl ;
    cout << q << " in HG00733" << endl ;
}

void Mixer::load_kmers(Track track) {
    auto c = Configuration::getInstance() ;
    int s = 0 ;
    bool found = false ;
    //cout << track.get_name() << endl ;
    for(auto sample: c->samples) {
        if (calls[s].find(track) != calls[s].end()) {
            bool correct = true ;
            //string g_1 = calls[s].find(track)->first.genotype ;
            //if (tracks[s].find(track) != tracks[s].end()) {
            //    string g_2 = tracks[s].find(track)->first.genotype ;
            //    if (g_1 == g_2 || (g_1.find("1") != string::npos && g_2.find("1") != string::npos)) {
            //        correct = true ;
            //    }
            //} else {
            //    if (g_1 == "0/0") {
            //        correct = true ;
            //    }
            //}
            if (correct) {
                //cout << "Genotype correct" << endl ;
                string path = c->workdir + "/" + sample + "/batch_" + track.get_name() + ".json" ;
                ifstream json_file(path) ;
                nlohmann::json kmers_json ;
                json_file >> kmers_json ;
                unordered_map<uint64_t, Kmer> _kmers ;
                for (nlohmann::json::iterator kmer = kmers_json.begin(); kmer != kmers_json.end(); kmer++) {
                    uint64_t k = encode_kmer(canonicalize(kmer.key()).c_str()) ;
                    Kmer _kmer = kmer.value().get<Kmer>() ;
                    _kmer.seq = k ; 
                    _kmer.count = 0 ;
                    _kmer.total = 0 ;
                    _kmers[k] = _kmer ;
                }
                if (_kmers.size() == 0) {
                    cout << track.get_name() << " " << sample << endl ;
                    cin.get() ;
                }
                json_file.close() ;
                if (!flags[track]) {
                    flags[track] = true ;
                    kmers[track].insert(_kmers.begin(), _kmers.end()) ;
                    //if (kmers[track].size() == 0) {
                    //    cout << "No kmers remainig for " << track.get_name() << endl ;
                    //}
                } else {
                    auto kmer = kmers[track].begin() ;
                    while (kmer != kmers[track].end()) {
                        if (_kmers.find(kmer->first) == _kmers.end()) {
                            kmer = kmers[track].erase(kmer) ;
                        } else {
                            kmer++ ;
                        }
                    }
                    //if (kmers[track].size() == 0) {
                    //    cout << "No kmers remainig for " << track.get_name() << endl ;
                    //}
                }
                //cout << kmers[track].size() << endl ;
            }
        } else {
            //cout << "Not found" << endl ;
        }
        s += 1 ;
    }
    //if (kmers[track].size() == 0) {
    //    cout << "No kmers remainig for " << track.get_name() << endl ;
    //}
}

void Mixer::export_kmers() {
    cout << "Exporting kmers.." << endl ;
    auto c = Configuration::getInstance() ;
    string path = c->workdir + "/Mix" ;
    vector<int> counters ;
    int num_batches = 100 ;
    vector<mutex> locks(num_batches) ;
    vector<ofstream> output_files ;
    for (int i = 0; i < num_batches; i++) {
        string p = path + "/kmers_batch_" + std::to_string(i) + ".json" ;
        output_files.emplace_back(ofstream {p}) ;
        output_files[i] << "{\n" ;
        counters.push_back(0) ;
    }
    int t = 0 ;
    int n = 0 ;
    for (auto track = kmers.begin(); track != kmers.end(); track++) {
        int m = track->second.size() ;
        n += m ;
        if (m != 0) {
            t += 1 ;
        }
    }
    cout << "Selected " << n << " kmers for " << t << " tracks.." << endl ;
    #pragma omp parallel for num_threads(c->threads)
    for (size_t bucket = 0; bucket < kmers.bucket_count(); bucket++) {
        for (auto track = kmers.begin(bucket); track != kmers.end(bucket); track++) {
            for (auto kmer = track->second.begin(); kmer != track->second.end(); kmer++) {
                int t = bucket % num_batches ;
                locks[t].lock() ;
                if (counters[t] != 0) {
                    output_files[t] << ",\n" ;
                }
                counters[t] += 1 ;
                output_files[t] << "\"" << decode_kmer(kmer->first) << "\":" ;
                output_files[t] << nlohmann::json(kmer->second).dump(4) ;
                locks[t].unlock() ;
            }
        }
    }
    for (int i = 0; i < num_batches; i++) {
        output_files[i] << "\n}\n" ;
    }
}
