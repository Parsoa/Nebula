#include <string>
#include <stdlib.h>

using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

#define MASK 6

// ============================================================================= \\
// ======================= Kmer Manipulation Helpers =========================== \\
// ============================================================================= \\

string reverse_complement(string s) {
    char rc[s.length()] ;
    int j = 0 ;
    for (int i = s.length() - 1 ; i >= 0 ; i--) {
        if (s[i] == 'A') {
            rc[j] = 'T' ;
        }
        else if (s[i] == 'T') {
            rc[j] = 'A' ;
        }
        else if (s[i] == 'C') {
            rc[j] = 'G' ;
        }
        else if (s[i] == 'G') {
            rc[j] = 'C' ;
        }
        else if (s[i] == 'N') {
            rc[j] == 'N' ;
        }
        j += 1 ;
    }
    return std::string(rc) ;
}

string canonicalize(string s) {
    string rc = reverse_complement(s) ;
    return s.compare(rc) <= 0 ? s : rc ;
}

uint64_t encode_kmer(const char* c) {
    uint64_t kmer = 0 ;
    for (uint64_t i = 0 ; i < 32 ; i++) {
        uint64_t d = (uint64_t)((toupper(c[i]) & MASK) >> 1) ;
        kmer += (d << (62 - (i * 2))) ;
    }
    return kmer ;
}

uint64_t encode_kmer(const string s) {
    return encode_kmer(s.c_str()) ;
}

uint64_t encode_substring(const char* c, int base, int l) {
    uint64_t mask = 0 ;
    for (uint64_t i = 0 ; i < l ; i++) {
        mask += ((uint64_t)((toupper(c[base + i]) & MASK) >> 1) << (62 - (i * 2))) ;
    }
    return mask ;
}

string decode_kmer(uint64_t kmer) {
    char* d = (char*) malloc(sizeof(char) * 33) ;
    d[32] = 0 ;
    uint64_t mask = 3 ;
    for (int i = 0; i < 32 ; i++) {
        uint64_t t = kmer & mask ;
        if (t == 0) {
            d[31 - i] = 'A' ;
        }
        else if (t == 1) {
            d[31 - i] = 'C' ;
        }
        else if (t == 2) {
            d[31 - i] = 'T' ;
        }
        else if (t == 3) {
            d[31 - i] = 'G' ;
        }
        kmer = kmer >> 2 ;
    }
    return std::string(d) ;
}

string decode_base(uint64_t base) {
    if (base == 0) {
        return "A" ;
    }
    else if (base == 1) {
        return "C" ;
    }
    else if (base == 2) {
        return "T" ;
    }
    else if (base == 3) {
        return "G" ;
    }
    return "N" ;
}

uint64_t MASK_Y = (uint64_t) 3 << 62 ;
uint64_t MASK_X = (uint64_t) 3 << 46 ;

bool is_subsequence(uint64_t x, uint64_t y) {
    int m = 0 ;
    int n = 0 ;
    int j = 46 ;
    uint64_t mask_y = MASK_Y ;
    uint64_t mask_x = MASK_X ;
    x = x >> 8 ;
    for (int i = 62; i >= 0 && j >= 0; i -= 2) {
        if ((x & mask_x) >> j == (y & mask_y) >> i) {
            j -= 2 ;
            mask_x = mask_x >> 2 ;
        } else {
            if (n < 2) {
                if ((x & (mask_x >> 2)) >> (j - 2) == (y & (mask_y >> 2)) >> (i - 2)) {
                    j -= 2 ;
                    mask_x = mask_x >> 2 ;
                }
            }
        }
        mask_y = mask_y >> 2 ; 
    }
    bool res = j == -2 ;
    return res ;
}

uint64_t encoded_reverse_complement(uint64_t kmer) {
    uint64_t rc = 0 ;
    uint64_t one = 1 ;
    uint64_t mask = (one << 62) + (one << 63) ;
    for (int i = 0; i < 32; i++) {
        uint64_t d = kmer & mask ;
        kmer = kmer << 2 ;
        d += (one << 63) ;
        rc += d ;
        if (i != 31) {
            rc = rc >> 2 ;
        }
    }
    return rc ;
}

