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
    void genotype_clusters() ;
    void genotype_cluster_lp(int, int) ;
    void genotype_cluster_likelihood(int batch_b, int batch_e) ;
    void export_genotypes() ;
    void filter_kmers() ;
    void remove_non_unique_kmers() ;
    
    double calc_kmer_genotype_likelihood(Kmer& kmer, std::vector<std::string> genotypes) ;
    void filter_kmer(Kmer& kmer, std::unordered_map<SimpleTrack, uint8_t>::iterator it, std::vector<std::string>* genotypes, double& max_likelihood, std::vector<std::string>* choice) ;

    std::unordered_map<uint64_t, Kmer> kmers ;
    std::unordered_map<uint64_t, Kmer> gc_kmers ;
    std::unordered_map<uint64_t, Kmer> depth_kmers ;

    double std ;
    double coverage ;
    std::vector<double> gc_coverage ;

    std::unordered_map<SimpleTrack, GenotypingTrack> genotyping_tracks ;
    std::unordered_map<SimpleTrack, int> cluster_index ;
    std::unordered_map<SimpleTrack, std::vector<uint64_t>> tracks ;
    std::vector<std::vector<SimpleTrack>> clusters ;

} ;

#endif
