#ifndef MIX_HPP
#define MIX_HPP

#include <vector>
#include <string>
#include <iostream>
#include <iterator>

#include "config.hpp"
#include "bed_utils.hpp"
#include "kmer_utils.hpp"
#include "chromosomes.hpp"

// mixes kmers from multiple preprocessing runs

class Mixer {

public:


private:

    std::unordered_map<uint64_t, Kmer> kmers ;
} ;

#endif
