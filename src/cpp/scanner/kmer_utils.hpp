#ifndef KMER_HPP
#define KMER_HPP

#include <string>
#include <stdlib.h>

// ============================================================================= \\
// ======================= Kmer Manipulation Helpers =========================== \\
// ============================================================================= \\

uint64_t encode_kmer(const char* c) ;
uint64_t encode_kmer(const string s) ;
uint64_t encoded_reverse_complement(uint64_t kmer) ;

string decode_kmer(uint64_t kmer) ;
string decode_base(uint64_t base) ;
string canonicalize(string s) ;
string reverse_complement(string s) ;

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

    KmerIterator(char* seq, int begin, int end): seq(seq), seq_begin(begin), seq_end(end), i(begin), done(false) {}

    explicit operator bool() const { return !done; }

    pointer operator->() const { return &kmer; }
    reference operator*() const { return kmer; }

    KmerIterator& operator++() ;
    KmerIterator operator++(int) ;

private:
    int i ;
    int seq_end ;
    int seq_begin ;
    bool done ;
    char* seq ;
    value_type kmer ;
};

#endif
