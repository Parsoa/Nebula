import os
import pwd
import sys
import argparse

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
    parser.add_argument("--bed", default = None)
    # the job to execute from this file
    parser.add_argument("--job")
    # path to a BED-like file containing a set of common SNVs to be considered when generating breakpoints
    parser.add_argument("--snp")
    # standard deviation to use for the normal distribution modeling kmers, separately calculated for each set of reads
    parser.add_argument("--std", type = int)
    # the seed to use for random number generation 
    parser.add_argument("--seed", type = int)
    # triggers the debug mode 
    parser.add_argument("--debug", action = 'store_true')
    # set of exons for coverage estimation 
    parser.add_argument("--exons")
    # the path to a previously generated khmer counttable without the .ct extension
    parser.add_argument("--khmer")
    # specifies that this counttable should return dummy values
    parser.add_argument("--dummy", action='store_true')
    # path to a FASTQ files, for jobs that need one as input
    parser.add_argument("--fastq", default = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref')
    # path to a FASTQ files, for jobs that need one as input
    parser.add_argument("--genes", default = '/share/hormozdiarilab/Codes/NebulousSerendipity/data/hgnc.txt')
    # the name of the genome being genotyped or whatver 
    parser.add_argument("--genome")
    # radius of the neighborhood considered for breakpoints 
    parser.add_argument("--radius", type = int, default = 50)
    # whether to generate random events during simulation or not
    parser.add_argument("--random", action = 'store_true')
    # whether to resume this job from reduce or not
    parser.add_argument("--reduce", action = 'store_true')
    # maximum number of cpu cores to use
    parser.add_argument("--threads", type = int, default = 48)
    # expected depth of coverage for the FASTQ file
    parser.add_argument("--coverage", type = int,default = 30)
    # the path to a jellyfish generated kmer count index
    parser.add_argument("--jellyfish", nargs = '*')
    # a reference genome assembly, used to extract sequences from a set of BED tracks etc
    parser.add_argument("--reference", default = 'hg19')
    # the outer insert size of the paired end reads 
    parser.add_argument("--insertsize", type = int, default = 1000)
    # indicates if this is part of a simulation
    parser.add_argument("--simulation", action = 'store_true')
    # the size of the reads in the FASTQ file 
    parser.add_argument("--readlength", type = int, default = 100)
    # whether the simulation should be heterozygous
    parser.add_argument("--heterozygous", action = 'store_true')
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
        seed = args.seed,
        debug = args.debug,
        exons = args.exons,
        khmer = args.khmer,
        ksize = 31,
        genome = args.genome,
        radius = args.radius,
        random = args.random,
        bed_file = args.bed,
        coverage =  1 * args.coverage if args.simulation else args.coverage,
        is_dummy = args.dummy,
        jellyfish = args.jellyfish, 
        fastq_file = os.path.abspath(args.fastq),
        simulation = args.simulation,
        genome_hg19 = genome_hg19,
        genome_hg38 = genome_hg38,
        max_threads = 1 if args.debug else args.threads,
        insert_size = args.insertsize,
        read_length = args.readlength,
        heterozygous = args.heterozygous,
        reference_genome = reference_genome,
        khmer_num_tables = khmer_num_tables,
        khmer_table_size = khmer_table_size,
        output_directory = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output')),
        resume_from_reduce = args.reduce,
        count_server_port = 6985,
    )
