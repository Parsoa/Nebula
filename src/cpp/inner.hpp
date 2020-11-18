#ifndef INN_HPP
#define INN_HPP

#include <vector>
#include <string>
#include <unordered_map>

#include "config.hpp"
#include "bed_utils.hpp"
#include "kmer_utils.hpp"

class InnerKmerExtractor {

public:

    InnerKmerExtractor(int) ;
    void run() ;
    std::unordered_map<uint64_t, Kmer> inner_kmers ;

private:
    
    std::unordered_map<uint64_t, Kmer> extract_deletion_inner_kmers(Track& track) ;
    std::unordered_map<uint64_t, Kmer> extract_insertion_inner_kmers(Track& track) ;
    void preprocess_inner_kmers() ;
    
    int threads ;
} ;

#endif
