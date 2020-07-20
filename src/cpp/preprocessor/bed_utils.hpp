#ifndef BED_HPP
#define BED_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <unordered_map>

#define SVTYPE_MISC         uint8_t(0) 
#define SVTYPE_DEL          uint8_t(1)
#define SVTYPE_INS          uint8_t(2)

struct Track {
    uint8_t chrom ;
    uint32_t begin ;
    uint32_t end ;
    uint8_t svtype ;
    uint32_t svlen ;
    std::string seq ;
} ;

uint8_t get_chromosome_index(std::string chrom) ;
std::string get_chromosome_name(uint8_t index) ;
int parse_svtype(std::string svtype) ;
Track parse_track_name(std::string name) ;
std::vector<Track> load_tracks_from_file(std::string path) ;

#endif
