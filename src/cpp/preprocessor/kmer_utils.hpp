#ifndef KMER_HPP
#define KMER_HPP

#include <vector>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <unordered_map>

#include "bed_utils.hpp"

#include "json.hpp"

#define KMER_TYPE_JUNCTION  uint8_t(0)
#define KMER_TYPE_INNER     uint8_t(1)
#define KMER_TYPE_DEPTH     uint8_t(2)
#define KMER_TYPE_GC        uint8_t(3)

#define LOCUS_TYPE_JUNCTION uint8_t(0)
#define LOCUS_TYPE_INNER    uint8_t(1)
#define LOCUS_TYPE_REF      uint8_t(2)

#define ITER_MODE_REF       uint8_t(0)
#define ITER_MODE_READ      uint8_t(1)

#define MASK 6

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

struct Locus {
    uint8_t chrom ;
    uint32_t position ;
    uint8_t type ; 
    uint64_t left ;
    uint64_t right ;
    uint16_t gc ;

    bool operator==(const Locus &o) const {
        return chrom == o.chrom && position == o.position && type == o.type ;
    }

    std::string get_name() const ;

    void parse_name(std::string name) ; 
} ;

namespace std {
    template <> struct hash<Locus> {
        std::size_t operator()(const Locus& l) const {
            uint64_t hash = 0 ;
            hash += l.position ;
            hash += (uint64_t) l.chrom << 32 ;
            hash += (uint64_t) l.type << 40 ;
            return std::hash<uint64_t>()(hash) ;
        }
    };
}

struct Kmer {
    uint8_t type ;
    uint32_t count ;
    uint32_t total ;
    uint32_t reference ;
    uint64_t seq ;
    std::unordered_map<SimpleTrack, uint8_t> tracks ;
    std::vector<Locus> loci ;
    std::vector<Locus> filtered_loci ;
    std::vector<Locus> junction_loci ;
    bool inverse ;
    uint16_t gc ;
    bool trend ;
    double weight ;

    Kmer(): inverse(false), count(0), total(0), reference(0), seq(0), type(KMER_TYPE_JUNCTION), weight(1.0) {}
    Kmer(uint64_t seq, uint64_t type): inverse(false), count(0), total(0), reference(0), seq(seq), type(type), weight(1.0) {}

    std::unordered_map<std::string, int> jsonify_tracks() const {
        std::unordered_map<std::string, int> _tracks ;
        for (auto track = tracks.begin(); track != tracks.end(); track++) {
            _tracks[track->first.get_name()] = track->second ;
        }
        return _tracks ;
    }
} ;

std::unordered_map<std::string, Locus> jsonify_loci(std::vector<Locus> loci) ;

void to_json(nlohmann::json& j, const Kmer& k) ;
void from_json(const nlohmann::json& j, Kmer& k) ;
void to_json(nlohmann::json& j, const Locus& l) ;
void from_json(const nlohmann::json& j, Locus& l) ;

// ============================================================================= \\
// ======================= Kmer Manipulation Helpers =========================== \\
// ============================================================================= \\

int calc_gc_content(const char* seq) ;

uint64_t encode_kmer(const char* c) ;
uint64_t encode_kmer(const std::string s) ;
uint64_t encoded_canonicalize(uint64_t kmer) ;
uint64_t encoded_reverse_complement(uint64_t kmer) ;

std::string decode_kmer(uint64_t kmer) ;
std::string decode_base(uint64_t base) ;
std::string canonicalize(std::string s) ;
std::string reverse_complement(std::string s) ;

bool is_subsequence(uint64_t x, uint64_t y) ;
bool is_canonical_subsequence(uint64_t x, uint64_t y) ;

struct SimpleKmer {
    uint64_t kmer ;
    uint64_t left ;
    uint64_t right ;
    uint16_t gc ;
};

class KmerIterator: public std::iterator<std::input_iterator_tag, SimpleKmer> {
public:
    // C++11 (explicit aliases)
    using iterator_category = std::input_iterator_tag ;
    using value_type = SimpleKmer ;
    using reference = value_type const& ;
    using pointer = value_type const* ;
    using difference_type = ptrdiff_t ;
    int position ;

    KmerIterator(const char* seq, int begin, int end, int position, int mode) ;

    explicit operator bool() const { return !done; }

    pointer operator->() const { return &kmer; }
    reference operator*() const { return kmer; }

    KmerIterator& operator++() ;
    KmerIterator operator++(int) ;

    void next() ;

    int l ;
    int i ; // current position in seq
    int seq_end ; // first kmer in seq
    int seq_begin ; // last kmer in seq
    int mode = ITER_MODE_REF ;
    bool done ;
    bool has_left = false ;
    bool has_right = false ;
    const char* seq ;
    value_type kmer ;
};

#endif
