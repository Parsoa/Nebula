#ifndef KMER_HPP
#define KMER_HPP

#include <string>
#include <stdlib.h>

// ============================================================================= \\
// ======================= Kmer Manipulation Helpers =========================== \\
// ============================================================================= \\

uint64_t encode_kmer(const char* c) ;
uint64_t encode_kmer(const std::string s) ;
uint64_t encoded_reverse_complement(uint64_t kmer) ;

std::string decode_kmer(uint64_t kmer) ;
std::string decode_base(uint64_t base) ;
std::string canonicalize(std::string s) ;
std::string reverse_complement(std::string s) ;

bool is_subsequence(uint64_t x, uint64_t y) ;
bool is_canonical_subsequence(uint64_t x, uint64_t y) ;

#endif
