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
                        snp,\
                        std,\
                        ksize,\
                        bed_file,\
                        coverage,\
                        is_dummy,\
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
            self.snp = snp
            self.std = std
            self.ksize = ksize
            self.bed_file = bed_file
            self.coverage = coverage
            self.is_dummy = is_dummy
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
    parser.add_argument("--bed", default = '/share/hormozdiarilab/Codes/NebulousSerendipity/data/CHM1_Lumpy.Del.100bp.DEL.bed')
    # path to a BED-like file containing a set of common SNVs to be considered when generating breakpoints
    parser.add_argument("--snp", default = '/share/hormozdiarilab/Codes/NebulousSerendipity/data/hg19.Common_SNPs.bed')
    # standard deviation to use for the normal distribution modeling kmers, separately calculated for each set of reads
    parser.add_argument("--std", type = int, default = 20)
    # specifies that this counttable should return dummy values
    parser.add_argument("--dummy", action='store_true')
    # path to a FASTQ files, for jobs that need one as input
    parser.add_argument("--fastq", default = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref')
    # maximum number of cpu cores to use
    parser.add_argument("--threads", type = int, default = 48)
    # expected depth of coverage for the FASTQ file
    parser.add_argument("--coverage", type = int, default = 30)
    # a reference genome assembly, used to extract sequences from a set of BED tracks etc
    parser.add_argument("--reference", default = 'hg19')
    # a FASTA/FASTQ/SAM/BAM file that should be used as the source for creating a counttable
    parser.add_argument("--counttable", default = '/share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio/HG00732.fq')
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
        snp = args.snp,
        std = args.std,
        ksize = 31,
        bed_file = os.path.abspath(args.bed),
        coverage = args.coverage,
        is_dummy = args.dummy,
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
 
