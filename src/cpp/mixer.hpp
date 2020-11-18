#ifndef MIX_HPP
#define MIX_HPP

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <unordered_map>

#include "bed_utils.hpp"
#include "kmer_utils.hpp"

// mixes kmers from multiple preprocessing runs

class Mixer {

public:

    void run() ;
    void load_tracks() ;
    void export_kmers() ;
    void load_kmers(Track track) ;

private:

    std::vector<std::unordered_map<Track, int>> calls ;
    std::vector<std::unordered_map<Track, int>> tracks ;
    std::unordered_map<Track, bool> flags ;
    std::unordered_map<Track, std::unordered_map<uint64_t, Kmer>> kmers ;

} ;

#endif
