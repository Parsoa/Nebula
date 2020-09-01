#ifndef CNT_HPP
#define CNT_HPP

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <unordered_map>

#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"

#include "config.hpp"
#include "bed_utils.hpp"
#include "kmer_utils.hpp"

#include "kseq.h"

#define COUNTER_MODE_BAM 0
#define COUNTER_MODE_FASTQ 1

KSEQ_INIT(gzFile, gzread)

class KmerCounter {

public:

    KmerCounter(int threads, std::unordered_map<uint64_t, Kmer>* gc_kmers, std::unordered_map<uint64_t, Kmer>* depth_kmers, std::unordered_map<uint64_t, Kmer>* genotyping_kmers): threads(threads), gc_kmers(gc_kmers), depth_kmers(depth_kmers), genotyping_kmers(genotyping_kmers) {}

    void run() ;
    void load_counts() ;

private:

    void load_kmers() ;
    void output_counts() ;
    void verify_counts() ;

    void process_read(const char* seq, int l, std::unordered_map<uint64_t, int>& _counts, std::unordered_map<uint64_t, int>& _totals) ;
    void process_reads(std::string) ;

    bool load_batch_bam(int p) ;
    void process_batch_bam(std::vector<bam1_t*> alignments, int thread, int p) ;

    bool load_batch_fastq(int p) ;
    void process_batch_fastq(std::vector<std::string> fastq_entries, int thread, int p) ;

    std::unordered_map<uint64_t, int> counts ;
    std::unordered_map<uint64_t, int> totals ;
    std::unordered_map<uint64_t, std::vector<uint64_t>*> masks ;

    std::vector<std::vector<std::unordered_map<uint64_t, int>>> partial_counts ;
    std::vector<std::vector<std::unordered_map<uint64_t, int>>> partial_totals ;

    gzFile fastq_file ;
    kseq_t* fastq_iterator ;
    std::vector<std::vector<std::vector<std::string>>> fastq_entries ;

    samFile *bam_file ;
    bam_hdr_t *bam_header ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;

    int mode ;
    int threads ;
    int batch_size ;

    std::unordered_map<uint64_t, Kmer>* gc_kmers ;
    std::unordered_map<uint64_t, Kmer>* depth_kmers ;
    std::unordered_map<uint64_t, Kmer>* genotyping_kmers ;
};

#endif
