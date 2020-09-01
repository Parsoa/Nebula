#include <string>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_map>

#include "mixer.hpp"
#include "config.hpp"
#include "genotyper.hpp"
#include "bed_utils.hpp"
#include "kmer_utils.hpp"
#include "chromosomes.hpp"
#include "preprocessor.hpp"

using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void load_tracks() {
    auto c = Configuration::getInstance() ;
    std::unordered_map<Track, int> tracks ;
    for (auto bed_file = c->bed.begin(); bed_file != c->bed.end(); bed_file++) {
        auto _tracks = load_tracks_from_file(*bed_file) ;
        for (auto track = _tracks.begin(); track != _tracks.end(); track++) {
            tracks[*track] = 1 ;
        }
    }
    for (auto track = tracks.begin(); track != tracks.end(); track++) {
        auto t = track->first ;
        SimpleTrack s{t.chrom, t.begin, t.end, t.svtype} ;
        c->tracks[s] = track->first ;
        bed_tracks.push_back(track->first) ;
    }
    cout << bed_tracks.size() << " tracks in total. " << endl ;
}

void create_workdir() {
    auto c = Configuration::getInstance() ;
    struct stat info ;
    if (stat(c->workdir.c_str(), &info) != 0) {
        cout << "Working directory does not exist. Aborting.." << endl ;
        exit(0) ;
    }
}

void test_print() {
    int threads = 24 ;
    std::mutex cout_mutex ;
    for (int t = 0; t < threads; t++) {
        cout << endl ;
    }
    cout << endl ;
    std::vector<uint64_t> counts(threads) ;
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < threads; i++) {
        while (true) {
            counts[i] += 1 ;
            if (counts[i] % (1000 + i) == 0) {
                cout_mutex.lock() ;
                for (int j = 0; j <= threads - i; j++) {
                    cout << "\x1b[A" ;
                }
                cout << "\r" ;
                cout << "Thread " << i << " printing.." ;
                for (int j = 0; j < threads - i; j++) {
                    cout << endl ;
                }
                cout_mutex.unlock() ;
            }
        }
    }
}

void unordered_map_test() {
    unordered_map<uint64_t, Kmer> kmers ;
    for (uint64_t i = 0; i < 10000000; i++) {
        kmers[i].seq = 0 ;
    }
    int n = 0 ;
    cout << kmers.bucket_count() << " buckets.. " << endl ;
    size_t m = 0 ;
    #pragma omp parallel for num_threads(48)
    for (size_t b = 0; b < kmers.bucket_count(); b++) {
        if (b > m) {
            m = b ;
        }
        for (auto kmer = kmers.begin(b); kmer != kmers.end(b); kmer++) {
            n += 1 ;
        }
        if (n % 10000 == 0) {
            cout << "processed " << n << " " << m << " " << kmers.bucket_count() << endl ;
        }
    }
}

int main(int argc, char** argv) {
    cout << "Nebula, ultra-efficient, mapping-free genotyper." << endl ;
    //test_print() ;
    //unordered_map_test() ;
    auto c = Configuration::getInstance() ;
    c->parse(argc - 1, argv + 1) ;
    create_workdir() ;
    load_tracks() ;
    if (strcmp(argv[1], "mix") == 0) {
        auto mixer = new Mixer() ;
        mixer->run() ;
    }
    if (strcmp(argv[1], "genotype") == 0) {
        auto genotyper = new Genotyper() ;
        genotyper->run() ;
    }
    if (strcmp(argv[1], "preprocess") == 0) {
        auto preprocessor = new Preprocessor() ;
        preprocessor->run() ;
    }
}

