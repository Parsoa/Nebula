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
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)

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
    # the job to execute from this file
    parser.add_argument("--job")
    # path to a BED-like file containing a set of common SNVs to be considered when generating breakpoints
    parser.add_argument("--snp")
    # standard deviation to use for the normal distribution modeling kmers, separately calculated for each set of reads
    parser.add_argument("--std", type = int)
    # specifies that this counttable should return dummy values
    parser.add_argument("--dummy", action='store_true')
    # path to a FASTQ files, for jobs that need one as input
    parser.add_argument("--fastq", default = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref')
        # path to a FASTQ files, for jobs that need one as input
    parser.add_argument("--genes", default = '/share/hormozdiarilab/Codes/NebulousSerendipity/data/hgnc.txt')
    # whether to resume this job from reduce or not
    parser.add_argument("--reduce", action = 'store_true')
    # maximum number of cpu cores to use
    parser.add_argument("--threads", type = int, default = 48)
    # expected depth of coverage for the FASTQ file
    parser.add_argument("--coverage", type = int)
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
    reference_genome = genome_hg19 if args.reference == 'hg19' else genome_hg38 if args.reference == 'hg38' else args.reference
    max_threads = args.threads
    khmer_table_size = 16e9
    khmer_num_tables = 4
    # set up configuration
    Configuration(
        job = args.job,
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
        resume_from_reduce = args.reduce,
        sample_count_server_port = 6985,
        reference_count_server_port = 8569
    )
    colorama.init()
 
