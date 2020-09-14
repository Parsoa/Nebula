#ifndef GNT_HPP
#define GNT_HPP

#include <vector>
#include <string>
#include <iostream>
#include <iterator>

#include "config.hpp"
#include "kmer_utils.hpp"

class Genotyper {

public:

    void run() ;

private:

    void genotype() ;
    void load_kmers() ;
    void load_counts() ;
    void count_kmers() ;
    void cluster_kmers() ;
    void estimate_coverage() ;
    void adjust_gc_coverage() ;
    void solve_lp() ;
    void genotype_clusters(int, int) ;
    void export_genotypes() ;
    void filter_kmers() ;

    std::unordered_map<uint64_t, Kmer> kmers ;
    std::unordered_map<uint64_t, Kmer> gc_kmers ;
    std::unordered_map<uint64_t, Kmer> depth_kmers ;

    double std ;
    double mean ;
    double coverage ;
    std::vector<double> gc_coverage ;

    std::unordered_map<SimpleTrack, GenotypingTrack> genotyping_tracks ;

    std::unordered_map<SimpleTrack, std::vector<uint64_t>> tracks ;
    std::unordered_map<SimpleTrack, int> cluster_index ;
    std::vector<std::vector<SimpleTrack>> clusters ;

} ;

#endif
