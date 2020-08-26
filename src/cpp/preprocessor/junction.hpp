#ifndef JNT_HPP
#define JNT_HPP

#include <vector>
#include <string>
#include <unordered_map>

#include "config.hpp"
#include "bed_utils.hpp"
#include "kmer_utils.hpp"

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"

class JunctionKmerExtractor {

public:

    JunctionKmerExtractor(int) ;
    
    void run() ;
    
    std::unordered_map<uint64_t, Kmer> junction_kmers ;

private:
    
    void load_tracks() ;
    std::unordered_map<uint64_t, Kmer> extract_junction_kmers(Track track, int t, int j) ;
    void preprocess_junction_kmers() ;
    
    int overlap(std::pair<int, int> a, std::pair<int, int> b) ;
    bool is_clipped(std::vector<std::pair<int, int>> clips, std::pair<int, int> kmer) ;

    int threads ;
    std::vector<std::vector<samFile*>> files ;
    std::vector<std::vector<bam_hdr_t*>> headers ;
    std::vector<std::vector<hts_idx_t*>> indices ;
    std::vector<std::vector<hts_itr_t*>> iterators ;

    std::vector<std::unordered_map<Track, int>> tracks ;
} ;

#endif
