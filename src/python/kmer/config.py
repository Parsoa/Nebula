import os
import pwd
import sys
import argparse

import colorama

print('importing config.py')
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Configuration:

    kmer_cache_size = 10000

    class __impl:
        def __init__(self,
                        ksize,\
                        bed_file,\
                        coverage,\
                        counttable,\
                        fastq_file,\
                        genome_hg19,\
                        genome_hg38,\
                        max_threads,\
                        reference_genome,\
                        khmer_num_tables,\
                        khmer_table_size,\
                        output_directory,\
                        sample_count_server_port,\
                        reference_count_server_port):
            self.ksize = ksize
            self.bed_file = bed_file
            self.coverage = coverage
            self.counttable = counttable
            self.fastq_file = fastq_file
            self.genome_hg19 = genome_hg19
            self.genome_hg38 = genome_hg38
            self.max_threads = max_threads
            self.reference_genome = reference_genome
            self.khmer_num_tables = khmer_num_tables
            self.khmer_table_size = khmer_table_size
            self.output_directory = output_directory
            self.count_server_port = sample_count_server_port
            self.sample_count_server_port = sample_count_server_port
            self.reference_count_server_port = reference_count_server_port

    __instance = None

    def __init__(self, **kwargs):
        if Configuration.__instance is None:
            Configuration.__instance = Configuration.__impl(**kwargs)

    def __getattr__(self, attr):
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):
        return setattr(self.__instance, attr, value)

# ============================================================================================================================ #
# Configuration
# ============================================================================================================================ #

def init():
    print('configuring...')
    configure(parse_args())
    print('configuration created')

def parse_args():
    parser = argparse.ArgumentParser()
    # path to a BED files, for jobs that need one as input
    parser.add_argument("--bed", default = 'CHM1_Lumpy.Del.100bp.DEL.bed')
    # path to a FASTQ files, for jobs that need one as input
    parser.add_argument("--fastq", default = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref')
    # maximum number of cpu cores to use
    parser.add_argument("--threads", type = int, default = 48)
    # expected depth of coverage for the FASTQ file
    parser.add_argument("--coverage", type = int, default = 30)
    # a reference genome assembly, used to extract sequences from a set of BED tracks etc
    parser.add_argument("--reference", default = 'hg38')
    # a FASTA/FASTQ/SAM/BAM file that should be used as the source for creating a counttable
    parser.add_argument("--counttable")
    args = parser.parse_args()
    #
    return args

def configure(args):
    genome_hg19 = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref'
    genome_hg38 = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa'
    max_threads = args.threads
    khmer_table_size = 16e9
    khmer_num_tables = 4
    #
    reference_genome = genome_hg19 if args.reference == 'hg19' else genome_hg38
    if not os.path.isfile(reference_genome):
        print("fatal error: couldn't find reference genome", args.reference, " aborting ...")
        exit()
    Configuration(
        ksize = 31,
        bed_file = os.path.abspath(args.bed),
        coverage = args.coverage,
        counttable = os.path.abspath(args.counttable),
        fastq_file = os.path.abspath(args.fastq),
        genome_hg19 = genome_hg19,
        genome_hg38 = genome_hg38,
        max_threads = max_threads, 
        reference_genome = reference_genome,
        khmer_num_tables = khmer_num_tables,
        khmer_table_size = khmer_table_size,
        output_directory = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output')),
        sample_count_server_port = 6985,
        reference_count_server_port = 8569
    )
    colorama.init()
 
