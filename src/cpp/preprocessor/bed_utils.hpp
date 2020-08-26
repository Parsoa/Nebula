#ifndef BED_HPP
#define BED_HPP

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#define SVTYPE_DEL          uint8_t(1)
#define SVTYPE_INS          uint8_t(2)
#define SVTYPE_MISC         uint8_t(3)

struct SimpleTrack {
    uint8_t chrom ;
    uint32_t begin ;
    uint32_t end ;
    uint8_t svtype ;

    bool operator==(const SimpleTrack &o) const {
        return chrom == o.chrom && begin == o.begin && end == o.end && svtype == o.svtype ;
    }

    bool operator!=(const SimpleTrack &o) const {
        return !(*this == o) ;
    }

    std::string get_name() const ;
    void parse_from_name(std::string name) ;
};

struct Track: SimpleTrack {
    uint32_t svlen ;
    std::string seq ;
    std::string genotype ;

    bool operator==(const Track &o) const {
        return chrom == o.chrom && begin == o.begin && end == o.end && svtype == o.svtype ;
    }
};

struct GenotypingTrack {
    double lp_value ;
};

namespace std {
    template <> struct hash<Track> {
        std::size_t operator()(const Track& t) const {
            uint64_t hash = 0 ;
            hash += t.begin ;
            hash += (uint64_t) t.chrom << 32 ;
            hash += (uint64_t) t.svtype << 40 ;
            return std::hash<uint64_t>()(hash) ;
        }
    };
    template <> struct hash<SimpleTrack> {
        std::size_t operator()(const SimpleTrack& t) const {
            uint64_t hash = 0 ;
            hash += t.begin ;
            hash += (uint64_t) t.chrom << 32 ;
            hash += (uint64_t) t.svtype << 40 ;
            return std::hash<uint64_t>()(hash) ;
        }
    };
}

extern std::vector<Track> bed_tracks ;

std::string get_svtype(uint8_t svtype) ;
uint8_t parse_svtype(std::string svtype) ;
std::string get_chromosome_name(uint8_t index) ;
uint8_t get_chromosome_index(std::string chrom) ;
std::vector<Track> load_tracks_from_file(std::string path) ;
std::unordered_map<Track, int> load_tracks_from_file_as_dict(std::string path) ;

#endif
