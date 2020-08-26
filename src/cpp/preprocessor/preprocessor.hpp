#ifndef PRE_HPP
#define PRE_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iostream>
#include <iterator>

#include "config.hpp"
#include "bed_utils.hpp"
#include "kmer_utils.hpp"
#include "chromosomes.hpp"

class Preprocessor {

public:

    void run() ;

private:

    void filter_kmers() ;
    void scan_reference(int threads) ;
    void dump_kmers(std::string path) ;
    std::unordered_map<uint64_t, Kmer> scan_chromosome(std::string chrom, int) ;
    void merge_genotyping_kmers(std::unordered_map<uint64_t, Kmer> inner_kmers, std::unordered_map<uint64_t, Kmer> junction_kmers) ;

    std::mutex cout_mutex ;
    std::unordered_map<Track, int> bed_tracks ;
    std::unordered_map<uint64_t, Kmer> kmers ;
} ;

#endif
