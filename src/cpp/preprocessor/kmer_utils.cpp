#include <assert.h>

#include "kmer_utils.hpp"

using namespace std ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

//unordered_map<string, Locus> jsonify_loci(vector<Locus> loci) {
//    unordered_map<string, Locus> _loci ;
//    for (auto locus = loci.begin(); locus != loci.end(); locus++) {
//        _loci[locus->get_name()] = *locus ;
//    }
//    return _loci ;
//}

string Locus::get_name() const {
    if (type == LOCUS_TYPE_INNER) {
        return "inside_" + get_chromosome_name(chrom) + "_" + std::to_string(position) ;
    }
    if (type == LOCUS_TYPE_JUNCTION) {
        return "junction_" + get_chromosome_name(chrom) + "_" + std::to_string(position) ;
    }
    return get_chromosome_name(chrom) + "_" + std::to_string(position) ;
} 

void Locus::parse_name(string name) {
    stringstream ss(name) ;
    vector<string> tokens ;
    string token ;
    while (getline(ss, token, '_')) {
        tokens.push_back(token) ;
    }
    type = tokens[0] == "junction" ? LOCUS_TYPE_JUNCTION : LOCUS_TYPE_INNER ;
    chrom = get_chromosome_index(tokens[1].substr(tokens[1].find("@") + 1)) ;
    position = std::stoi(tokens[2]) ;
}

void to_json(nlohmann::json& j, const Kmer& k) {
    j = nlohmann::json{{"count", k.count}, {"reference", k.reference}, {"loci", k.loci},\
        {"junction_loci", k.junction_loci}, {"filtered_loci", k.filtered_loci}, {"total", k.total}, \
        {"tracks", k.jsonify_tracks()}, {"inverse", k.inverse}, {"trend", k.trend}, {"gc", k.gc}, {"type", k.type}
    } ;
}

void to_json(nlohmann::json& j, const Locus& l) {
    nlohmann::json masks({}) ;
    if (l.left != 0 ) {
        masks[decode_kmer(l.left)] = 0 ;
    }
    if (l.right != 0) {
        masks[decode_kmer(l.right)] = 0 ;
    }
    j = nlohmann::json{{"chrom", get_chromosome_name(l.chrom)}, {"position", l.position}, {"masks", masks}, {"type", l.type}, {"gc", l.gc}} ;
}

void from_json(const nlohmann::json& j, Locus& l) {
    nlohmann::json masks = j.at("masks") ;
    if (masks.find("left") != masks.end()) {
        l.left = encode_kmer(masks.at("left")) ;
    } else {
        l.left = 0 ;
    }
    if (masks.find("right") != masks.end()) {
        l.right = encode_kmer(masks.at("right")) ;
    } else {
        l.right = 0 ;
    }
}

void from_json(const nlohmann::json& j, Kmer& k) {
    if (j.find("type") != j.end()) {
        k.type = j["type"] ;
    }
    if (j.find("trend") != j.end()) {
        k.trend = j["trend"] ;
    }
    if (j.find("inverse") != j.end()) {
        k.inverse = j["inverse"] ;
    }
    if (j.find("reference") != j.end()) {
        k.reference = j["reference"] ;
    }
    for (auto track = j.at("loci").begin(); track != j.at("loci").end(); track++) {
        Locus l ;
        k.loci.push_back(l) ;
    }
    for (auto track = j.at("tracks").begin(); track != j.at("tracks").end(); track++) {
        SimpleTrack t ;
        t.parse_from_name(track.key()) ;
        k.tracks[t] = track.value() ;
    }
    for (auto track = j.at("junction_loci").begin(); track != j.at("junction_loci").end(); track++) {
        Locus l ;
        k.junction_loci.push_back(l) ;
    }
    for (auto track = j.at("filtered_loci").begin(); track != j.at("filtered_loci").end(); track++) {
        Locus l ;
        k.filtered_loci.push_back(l) ;
    }
}

// ============================================================================= \\
// ======================= Kmer Manipulation Helpers =========================== \\
// ============================================================================= \\

KmerIterator& KmerIterator::operator++() {
    assert(!done) ;
    if (i >= seq_end) {
        done = true ;
        return *this ;
    } else {
        next() ;
        return *this ;
    }
}
    
KmerIterator::KmerIterator(const char* seq, int begin, int end, int position, int mode): seq(seq), seq_begin(begin), seq_end(end), i(begin), position(position), done(false), mode(mode) {
    l = strlen(seq) ;
    if (l < 32 || end < begin) {
        done = true ;
    } else {
        if (mode == ITER_MODE_READ) {
            for (int j = 0; j < l; j++) {
                if (seq[j] == 'N') {
                    done = true ;
                }
            }
        }
        kmer.left = 0 ;
        kmer.right = 0 ;
        if (i >= 32) {
            has_left = true ;
            kmer.left = encode_kmer(seq + i - 32) ;
        }
        kmer.kmer = encode_kmer(seq + i) ;
        if (i + 32 + 32 <= strlen(seq)) {
            has_right = true ;
            kmer.right = encode_kmer(seq + i + 32) ;
        }
        if (i - 234 >= 0 && i - 234 + 500 <= l) {
            kmer.gc = calc_gc_content(seq + i - 234) ;
        }
        i++ ;
        position++ ;
    }
}

void KmerIterator::next() {
    kmer.kmer = kmer.kmer << 2 ;
    kmer.kmer += (seq[i + 31] & MASK) >> 1 ;
    // right mask
    if (i + 32 + 32 <= l) {
        kmer.right = kmer.right << 2 ;
        kmer.right += (seq[i + 32 + 31] & MASK) >> 1 ;
    } else {
        has_right = false ;
        kmer.right = 0 ;
    }
    // left mask
    if (i == 32) {
        has_left = true ;
        kmer.left = encode_kmer(seq + i - 32) ;
    } else if (i > 32) {
        kmer.left = kmer.left << 2 ;
        kmer.left += (seq[i - 1] & MASK) >> 1 ;
    } else {
        kmer.left = 0 ;
    }
    // GC content
    if (i - (234 + 1) >= 0) {
        if (seq[i - (234 + 1)] == 'C' || seq[i - (234 + 1)] == 'G') {
            kmer.gc -= 1 ;
        }
    }
    if (i + (234 - 1 + 32) < l) {
        if (seq[i + (234 - 1 + 32)] == 'C' || seq[i + (234 - 1 + 32)] == 'G') {
            kmer.gc += 1 ;
        }
    }
    i++ ;
    position++ ;
}

KmerIterator KmerIterator::operator++(int) {
    KmerIterator const tmp(*this) ;
    ++*this ;
    return tmp ;
}

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

// Checks is x is a subsequence of y
bool is_subsequence(uint64_t x, uint64_t y) {
    int m = 0 ;
    int n = 0 ;
    // ignore right 4 bases, 56 bits effective
    x = x >> 8 ;
    // ignore left 4 bases, 48 bases effective
    int j = 46 ;
    uint64_t mask_x = (uint64_t) 3 << 46 ;
    // all bases considered
    uint64_t mask_y = (uint64_t) 3 << 62 ;
    for (int i = 62; i >= 0 && j >= 0; i -= 2) {
        if ((x & mask_x) >> j == (y & mask_y) >> i) {
            j -= 2 ;
            mask_x = mask_x >> 2 ;
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

uint64_t encoded_canonicalize(uint64_t kmer) {
    uint64_t rc = encoded_reverse_complement(kmer) ;
    if (kmer <= rc) {
        return kmer ;
    }
    return rc ;
}

bool is_canonical_subsequence(uint64_t x, uint64_t y) {
    return is_subsequence(x, y) || is_subsequence(encoded_reverse_complement(x), y) ;
}

int calc_gc_content(const char* seq) {
    int gc = 0 ;
    for (int i = 0; i < 500; i++) {
        if (seq[i] == 'C' || seq[i] == 'G') {
            gc += 1 ;
        }
    }
    return gc / 5 ;
}

