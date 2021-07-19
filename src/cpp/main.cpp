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
        cout << "Working directory does not exist. creating.." << endl ;
        int check = mkdir(c->workdir.c_str(), 0777) ;
        if (check != 0) {
            cerr << "Failed to create output directory, aborting.." << endl ;
            exit(check) ;
        }
    }
}

int main(int argc, char** argv) {
    cerr << "Nebula, ultra-efficient, mapping-free genotyper." << endl ;
    if (argc == 1) { 
        cerr << "Basic usage:" << endl ;
        cerr << "\t - To extract kmers (preprocessing):" << endl ;
        cerr << "\t\t nebula preprocess --bed /path/to/genotypes.bed --bam /path/to/bam_file.bed --wokdir <path> --reference /path/to/reference/FASTA/file --thread <integer>" << endl; 
        cerr << "\t - To genotype a sample:" << endl ;
        cerr << "\t\t nebula genotype [--bed /path_to_genotypes.bed] --bam /path/to/bam_file.bed --workdir <path> --kmers /path/to/kmers --depth_kmers depth_kmers.json --gc_kmers gc_kmers.json [--select] [--unique]" << endl ;
        cerr << "\t - To mix kmers from a preprocessing run:" << endl ;
        cerr << "\t\t nebula mix --bed /path_to_genotypes_1.bed,/path_to_genotypes_2.bed --samples sample_1,sample_2 --workdir /output" << endl ;
        cerr << "\t This only details the list of options for each subcommand. See https://github.com/Parsoa/Nebula for detailed instructions." << endl ;
        exit(0) ;
    }
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

