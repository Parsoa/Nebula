#include <omp.h>
#include <ctime>
#include <string>
#include <thread>
#include <locale>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <pthread.h>
#include <unordered_map>

#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <hts_internal.h>
#include "htslib/hts_endian.h"

#include "htslib/hfile.h"

#include "kseq.h"
#include "json.hpp"

KSEQ_INIT(gzFile, gzread)
using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

#define MASK 6

std::unordered_map<uint64_t, int> counts ;
std::unordered_map<uint64_t, int> totals ;
std::unordered_map<uint64_t, std::vector<uint64_t>*> masks ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void __bam_cigar2rqlens(int n_cigar, const uint32_t *cigar, int *rlen, int *qlen)
{
    int k;
    *rlen = *qlen = 0;
    for (k = 0; k < n_cigar; ++k) {
        int type = bam_cigar_type(bam_cigar_op(cigar[k]));
        int len = bam_cigar_oplen(cigar[k]);
        if (type & 1) *qlen += len;
        if (type & 2) *rlen += len;
    }
}

static void __swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

static int __do_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
    //new_data = (uint8_t*) realloc(b->data, new_m_data);
    new_data = (uint8_t*) malloc(sizeof(uint8_t) * new_m_data);
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}

static inline int __realloc_bam_data(bam1_t *b, size_t desired)
{
    //if (desired <= b->m_data) return 0;
    return __do_realloc_bam_data(b, desired);
}

static inline int __possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    uint32_t new_len = b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < b->l_data) {
        errno = ENOMEM;
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return __do_realloc_bam_data(b, new_len);
}

static int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
{
    bam1_core_t *c = &b->core;
    uint32_t cigar_st, n_cigar4, CG_st, CG_en, ori_len = b->l_data, *cigar0, CG_len, fake_bytes;
    uint8_t *CG;

    // test where there is a real CIGAR in the CG tag to move
    if (c->n_cigar == 0 || c->tid < 0 || c->pos < 0) return 0;
    cigar0 = bam_get_cigar(b);
    if (bam_cigar_op(cigar0[0]) != BAM_CSOFT_CLIP || bam_cigar_oplen(cigar0[0]) != c->l_qseq) return 0;
    fake_bytes = c->n_cigar * 4;
    if ((CG = bam_aux_get(b, "CG")) == 0) return 0; // no CG tag
    if (CG[0] != 'B' || CG[1] != 'I') return 0; // not of type B,I
    CG_len = le_to_u32(CG + 2);
    if (CG_len < c->n_cigar || CG_len >= 1U<<29) return 0; // don't move if the real CIGAR length is shorter than the fake cigar length

    // move from the CG tag to the right position
    cigar_st = (uint8_t*)cigar0 - b->data;
    c->n_cigar = CG_len;
    n_cigar4 = c->n_cigar * 4;
    CG_st = CG - b->data - 2;
    CG_en = CG_st + 8 + n_cigar4;
    if (__possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
    b->l_data = b->l_data - fake_bytes + n_cigar4; // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes, ori_len - (cigar_st + fake_bytes)); // insert c->n_cigar-fake_bytes empty space to make room
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4); // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
    if (ori_len > CG_en) // move data after the CG tag
        memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
    if (recal_bin)
        b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)), 14, 5);
    if (give_warning)
        hts_log_error("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
    return 1;
}


int __bam_read1(BGZF *fp, bam1_t *b)
{
    //bam1_core_t *c = &b->core;
    bam1_core_t d ;
    bam1_core_t *c = &d ; //new bam1_core_t() ;
    b = new bam1_t() ;

    int32_t block_len, ret, i;
    uint32_t x[8], new_l_data;
    if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
        if (ret == 0) return -1; // normal end-of-file
        else return -2; // truncated
    }
    if (fp->is_be)
        ed_swap_4p(&block_len);
    if (block_len < 32) return -4;  // block_len includes core data
    if (bgzf_read(fp, x, 32) != 32) return -3;
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
    }
    c->tid = x[0]; c->pos = x[1];
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
    if ((uint32_t) c->l_qname + c->l_extranul > 255) // l_qname would overflow
        return -4;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];

    new_l_data = block_len - 32 + c->l_extranul;
    if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) return -4;
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
        return -4;
    //TODO reallocate for b
    if (__realloc_bam_data(b, new_l_data) < 0) return -4;
    b->l_data = new_l_data;

    if (bgzf_read(fp, b->data, c->l_qname) != c->l_qname) return -4;
    for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname+i] = '\0';
    c->l_qname += c->l_extranul;
    if (b->l_data < c->l_qname ||
        bgzf_read(fp, b->data + c->l_qname, b->l_data - c->l_qname) != b->l_data - c->l_qname)
        return -4;
    if (fp->is_be) __swap_data(c, b->l_data, b->data, 0);
    if (bam_tag2cigar(b, 0, 0) < 0)
        return -4;
    
    b->core = d ;
    if (c->n_cigar > 0) { // recompute "bin" and check CIGAR-qlen consistency
        int rlen, qlen;
        __bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
        if ((b->core.flag & BAM_FUNMAP)) rlen=1;
        b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
        // Sanity check for broken CIGAR alignments
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
            hts_log_error("CIGAR and query sequence lengths differ for %s",
                    bam_get_qname(b));
            return -4;
        }
    }
    return 4 + block_len;
}

int __sam_read1(htsFile *fp, bam_hdr_t *h, bam1_t *b) {
    int r = __bam_read1(fp->fp.bgzf, b);
    if (h && r >= 0) {
        if (b->core.tid  >= h->n_targets || b->core.tid  < -1 ||
            b->core.mtid >= h->n_targets || b->core.mtid < -1) {
            return -3;
        }
    }
    return r;
}

// ============================================================================= \\
// ============================================================================= \\
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
        uint64_t d = (uint64_t)((c[i] & MASK) >> 1) ;
        kmer += (d << (62 - (i * 2))) ;
    }
    return kmer ;
}

uint64_t encode_substring(const char* c, int base, int l) {
    uint64_t mask = 0 ;
    for (uint64_t i = 0 ; i < l ; i++) {
        mask += ((uint64_t)((c[base + i] & MASK) >> 1) << (62 - (i * 2))) ;
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

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void process_read(const char* seq, int l, unordered_map<uint64_t, int>& _counts) {
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    //cout << seq << endl ;
    for (int i = 0 ; i <= l - 32 ; i++) {
        if (i == 0) {
            k = encode_kmer(seq) ;
            left = encode_substring(seq, 32, 32) ;
            right = k ;
        } else {
            k = k << 2 ;
            k += (seq[i + 31] & MASK) >> 1 ;
        }
        if (i + 32 + 31 < l) {
            left = left << 2 ;
            left += (seq[i + 32 + 31] & MASK) >> 1 ;
        }
        if (i > 32) {
            right = right << 2 ;
            right += (seq[i - 1] & MASK) >> 1 ;
        }
        if (counts.find(k) != counts.end()) {
            totals[k] += 1 ;
            if (masks.find(k) == masks.end()) {
                if (_counts.find(k) == _counts.end()) {
                    _counts[k] = 0 ;
                }
                _counts[k] += 1 ;
            } else {
                std::vector<uint64_t>* m = masks[k] ;
                for (std::vector<uint64_t>::iterator it = m->begin(); it != m->end(); it++) {
                    if (i + 32 + 31 < l) {
                        if (is_subsequence(*it, left)) {
                            if (_counts.find(k) == _counts.end()) {
                                _counts[k] = 0 ;
                            }
                            _counts[k] += 1 ;
                            break ;
                        }
                    }
                    if (i >= 32) {
                        if (is_subsequence(*it, right)) {
                            if (_counts.find(k) == _counts.end()) {
                                _counts[k] = 0 ;
                            }
                            _counts[k] += 1 ;
                            break ;
                        }
                    }
                }
            }
        }
    }
}

// ============================================================================= \\
// ================================= BAM Files ================================= \\
// ============================================================================= \\

samFile *bam_file ;
bam_hdr_t *bam_header ;
bam1_t *bam_iterator ;

vector<vector<vector<bam1_t*>>> bam_entries ;

vector<samFile*> bam_files ;
vector<bam_hdr_t*> bam_headers ;
vector<bam1_t*> bam_iterators ;

bool load_batch_bam(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        bam_entries[p][i].clear() ;
    }
    int n = 0 ;
    while (sam_read1(bam_file, bam_header, bam_iterator) > 0) {
        //cout << bam_iterator << endl ;
        //uint8_t* data = (uint8_t*) malloc(sizeof(uint8_t) * bam_iterator->l_data) ;
        //memcpy(data, bam_iterator->data, bam_iterator->l_data) ;
        //bam1_t bam_entry ;
        //bam_entry.core = bam_iterator->core ;
        //bam_entry.id = bam_iterator->id ;
        //bam_entry.data = data ;
        //bam_entry.l_data = bam_iterator->l_data ;
        //bam_entry.m_data = bam_iterator->m_data ;
        //= {bam_iterator->core, bam_iterator->id, data, bam_iterator->l_data, bam_iterator->m_data} ;
        bam_entries[p][n % threads].push_back(bam_iterator) ;
        n += 1 ;
        if (n == batch_size) {
            return true ;
        }
        bam_iterator = bam_init1() ; 
    }
    if(bam_entries[p].empty()) {
        return false ;
    }
    return true ;
}

unordered_map<uint64_t, int> process_batch_bam(vector<bam1_t*> alignments) {
    char* seq ; //= (char*) malloc(200) ;
    uint32_t len = 0 ;
    unordered_map<uint64_t, int> _counts ;
    for (const auto alignment : alignments) {
        uint32_t l = alignment->core.l_qseq ; //length of the read
        if (l > len) {
            if (len > 0) {
                free(seq) ;
            }
            len = l ;
            seq = (char*) malloc(l + 1) ;
        }
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < len; i++){
            seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        seq[l] = '\0' ; // null terminate
        process_read(seq, l, _counts) ;
    }
    
    return _counts ;
}

//unordered_map<uint64_t, int> process_batch_bam(vector<bam1_t> alignments, int index, int threads, int batch_size) {
//    char* seq ;
//    uint32_t len = 0 ;
//    unordered_map<uint64_t, int> _counts ;
//    int n = 0 ;
//    int m = 0 ;
//    while (__sam_read1(bam_files[index], bam_headers[index], bam_iterators[index]) > 0) {
//        if (n % threads == index) {
//            //cout << n << endl ;
//            //cout << alignment << endl ;
//            uint32_t l = bam_iterators[index]->core.l_qseq ; //length of the read
//            if (l > len) {
//                if (len > 0) {
//                    free(seq) ;
//                }
//                len = l ;
//                seq = (char*) malloc(l + 1) ;
//            }
//            uint8_t *q = bam_get_seq(bam_iterators[index]) ; //quality string
//            for (int i = 0; i < len; i++){
//                seq[i] = seq_nt16_str[bam_seqi(q, i)] ; //gets nucleotide id and converts them into IUPAC id.
//            }
//            seq[l] = '\0' ; // null terminate
//            process_read(seq, l + 1, _counts) ;
//            //free(alignment.data) ;
//            m += 1 ;
//            if (m == batch_size / threads) {
//                break ;
//            }
//        }
//        n += 1 ;
//    }
//    return _counts ;
//}

// ============================================================================= \\
// ================================ FASTQ Files ================================ \\
// ============================================================================= \\

gzFile fastq_file ;
kseq_t* fastq_iterator ;
vector<vector<vector<string>>> fastq_entries ;

bool load_batch_fastq(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        fastq_entries[p][i].clear() ;
    }
    int l = 0 ;
    int n = 0 ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        fastq_entries[p][n % threads].push_back(fastq_iterator->seq.s) ;
        n += 1 ;
        if (n == batch_size) {
            return true ;
        }
    }
    if(fastq_entries[p].empty()) {
        return false ;
    }
    return true ;
}

unordered_map<uint64_t, int> process_batch_fastq(vector<string> fastq_entries) {
    int l ;
    const char* seq ;
    unordered_map<uint64_t, int> _counts ;
    for (const auto fastq_entry : fastq_entries) {
        seq = fastq_entry.c_str() ; 
        l = strlen(seq) ;
        process_read(seq, l, _counts) ;
    }
    return _counts ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

#define BAM 0
#define FASTQ 1

int process_reads(int threads, int mode, string input_file) {
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(threads)) ;
        fastq_entries.push_back(vector<vector<string>>(threads)) ;
    }
    for (int i = 0; i < threads; i++) {
        bam_files.push_back(hts_open(input_file.c_str(), "r")) ;
        bam_headers.push_back(sam_hdr_read(bam_files[i])) ;
        bam_iterators.push_back(bam_init1()) ;
    }
    int p = 0 ;
    int batch_size = 2000000 ;
    vector<unordered_map<uint64_t, int>> batches(threads) ;
    if (mode == BAM) {
        load_batch_bam(threads, batch_size, p) ;
    } else {
        load_batch_fastq(threads, batch_size, p) ;
    }
    cout << "Loaded batch" << endl ;
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;
    while (true) {
        uint64_t v = u ;
        for (int i = 0 ; i < threads ; i++) {
            if (mode == BAM) {
                u += bam_entries[p][i].size() ;
            } else {
                u += fastq_entries[p][i].size() ;
            }
        }
        if (v == u) {
            break ;
        }
        #pragma omp parallel for num_threads(threads + 1)
        for(int i = 0; i <= threads; i++) {
            if (i == 0) {
                if (mode == BAM) {
                    load_batch_bam(threads, batch_size, (p + 1) % 2, b) ;
                } else {
                    load_batch_fastq(threads, batch_size, (p + 1) % 2) ;
                }
            } else {
                unordered_map<uint64_t, int> _counts ;
                if (mode == BAM) {
                    //_counts = process_batch_bam(bam_entries[p][i - 1], i - 1, threads,  batch_size) ;
                    _counts = process_batch_bam(bam_entries[p][i - 1]) ;
                } else {
                    _counts = process_batch_fastq(fastq_entries[p][i - 1]) ;
                }
                batches[i - 1] = _counts ;
            }
        }
        p += 1 ;
        p %= 2 ;
        for (const auto &batch : batches) {
            for (const auto kmer : batch) {
                counts[kmer.first] += kmer.second ;
            }
        }
        b += 1 ;
        time_t s ;
        time(&s) ;
        cerr << "Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << u << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << "\r" ;
    }
    return u ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void output_counts(string path) {
    nlohmann::json payload ;
    cout << "Dumping kmer counts..." << endl ;
    string p = path + "/counts.json" ;
    cout << path << endl ;
    std::ofstream o(p);
    for (const auto kmer : counts) {
        o << decode_kmer(kmer.first) << ":" << kmer.second << ":" << totals[kmer.first] << "\n" ;
    }
}

int load_kmers(string path) {
    cout << "Loading kmers from " << path << ".." << endl ;
    ifstream json_file(path + "/pre_inner_kmers.json") ;
    nlohmann::json kmers_json ;
    json_file >> kmers_json ;
    for (nlohmann::json::iterator it = kmers_json.begin(); it != kmers_json.end(); ++it) {
        auto kmer = it.value() ;
        uint64_t k = encode_kmer(std::string(it.key()).c_str()) ;
        uint64_t rc_k = encode_kmer(reverse_complement(std::string(it.key())).c_str()) ;
        counts[k] = 0 ;
        counts[rc_k] = 0 ;
        totals[k] = 0 ;
        totals[rc_k] = 0 ;
        if (kmer["loci"].size()) {
            // keep this as a pointer for memory optimization, avoids keeping two copies
            std::vector<uint64_t> *m = new std::vector<uint64_t> ;
            for (nlohmann::json::iterator locus = kmer["loci"].begin(); locus != kmer["loci"].end(); ++locus) {
                for (nlohmann::json::iterator mask = locus.value()["masks"].begin(); mask != locus.value()["masks"].end(); ++mask) {
                    m->push_back(encode_kmer(mask.key().c_str())) ;
                    m->push_back(encode_kmer(reverse_complement(mask.key()).c_str())) ;
                }
            }
            if (m->size()) {
                masks[k] = m ;
                masks[rc_k] = m ;
            } 
        }
    }
    cout << "Loaded " << counts.size() / 2 << " kmers." << endl ;
    return 0 ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int main(int argc, char** argv) {
    string path(argv[1]) ;
    string input_file(argv[2]) ;
    int threads = std::stoi(string(argv[3]), nullptr, 10) ;
    time_t t ;
    time(&t) ;
    load_kmers(path) ;
    int u = 0;
    if (input_file.compare(input_file.size() - 4, 4, ".bam") == 0) {
        cout << "Input is BAM." << endl ;
        bam_file = hts_open(input_file.c_str(), "r") ;
        bam_header = sam_hdr_read(bam_file) ; //read header
        bam_iterator = bam_init1() ; //initialize an alignment
        u = process_reads(threads, BAM, input_file) ;
    } else {
        cout << "Input is FASTQ." << endl ;
        fastq_file = gzopen(input_file.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
        u = process_reads(threads, FASTQ, input_file) ;
    }
    time_t s ;
    time(&s) ;
    auto d = s - t ;
    cout << "\33[2K" << "Processed " << u << " reads in " << (s - t) << " seconds. Averaging " << u / d << " reads per second." << endl ;
    output_counts(path) ;
    cout << "Returning to CgcCounterJob.." << endl ;
}
