import os
import pwd
import sys
import argparse

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Configuration:

    class __impl:
        def __init__(self, args):
            for arg in args:
                for attr, value in arg.__dict__.items():
                    setattr(self, attr, value)
            setattr(self, 'hsize', self.ksize / 2)
            if hasattr(self, 'simulation'):
                if self.simulation:
                    self.simulation = int(self.simulation[:-1])
            if self.debug:
                self.threads = 1
            self.log_level = 1

    __instance = None
    __cmd_options = None

    @staticmethod
    def update(d):
        for key in d:
            setattr(Configuration.__instance, key, d[key])

    def __init__(self, *args):
        if Configuration.__instance is None:
            Configuration.__cmd_options = args
            Configuration.__instance = Configuration.__impl(args)

    def __getattr__(self, attr):
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):
        return setattr(self.__instance, attr, value)

# ============================================================================================================================ #
# Configuration
# ============================================================================================================================ #

def init():
    Configuration(parse_args())

def parse_args():
    parser = argparse.ArgumentParser(add_help = False)
    ################## General arguments
    # path to a BED files, for jobs that need one as input
    parser.add_argument("--bed", nargs = '+', default = None)
    # path to a VCF files, for jobs that need one as input
    parser.add_argument("--vcf", nargs = '+', default = None)
    # path to a BAM files, for jobs that need one as input
    parser.add_argument("--bam", default = None)
    # whether we are running on CGC or not 
    parser.add_argument("--cgc", action = 'store_true')
    # gap size, should be odd
    parser.add_argument("--gap", default = None, type = int) #TODO: depracate
    # the job to execute from this file
    parser.add_argument("--job")
    # whether we are running on rum or not
    parser.add_argument("--rum", action = 'store_true')
    # std of the kmer normal distribution
    parser.add_argument("--std", default = 15, type = int)
    # the seed to use for random number generation
    parser.add_argument("--seed", type = int, default = 0)
    # triggers the debug mode
    parser.add_argument("--debug", action = 'store_true')
    # set of exons for the current reference
    parser.add_argument("--exons", default = '/share/hormozdiarilab/Codes/NebulousSerendipity/data/Exons/hg38.exons.filtered.bed')
    # path to a FASTQ files, for jobs that need one as input
    parser.add_argument("--fastq", nargs = '*', default = None)
    # generic flag for passing input arguments
    parser.add_argument("--input")
    # a JSON file containing the list of kmers to count
    parser.add_argument("--kmers", nargs = '+')
    # length of the kmers 
    parser.add_argument("--ksize", default = 32, type = int)
    # assembled contigs for each structural variation
    parser.add_argument("--contigs", default = None)
    # the name of the genome being genotyped or whatver 
    parser.add_argument("--genome")
    # a JSON file containing the list of tracks and their kmers
    parser.add_argument("--tracks")
    # whether to resume this job from reduce or not
    parser.add_argument("--reduce", action = 'store_true')
    # which solver to use 
    parser.add_argument("--solver", default = 'coin')
    # maximum number of cpu cores to use
    parser.add_argument("--threads", type = int, default = 48)
    # alternate directory for previous job
    parser.add_argument("--previous", default = None)
    # expected depth of coverage for the FASTQ file
    parser.add_argument("--coverage", type = float, default = 50)
    # the path to a jellyfish generated kmer count index
    parser.add_argument("--jellyfish")
    # a reference genome assembly, used to extract sequences from a set of BED tracks etc
    parser.add_argument("--reference", default = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref')
    # tmp file directory 
    parser.add_argument("--workdir", default = '/share/hormozdiarilab/Codes/NebulousSerendipity/output')
    # path to a set of kmers for gc content estimation
    parser.add_argument("--gckmers", default = None)
    # if selecting kmers for genotyping 
    parser.add_argument("--select", action = 'store_true')
    ## Simulation options
    # whether to generate random events during simulation or not
    parser.add_argument("--random", action = 'store_true')
    # whether to do a diploid simulation or two haploid ones 
    parser.add_argument("--diploid", action = 'store_true')
    # indicates if this is part of a simulation
    parser.add_argument("--simulation", default = None)
    # if we only want to simulate a certain chromosome
    parser.add_argument("--chromosomes", default = None, nargs = '*')
    # the size of the reads in the fastq file 
    parser.add_argument("--readlength", type = int, default = 100)
    # description of this simulation 
    parser.add_argument("--description")
    # rate of SNPs in simulation 
    parser.add_argument("--mutation_rate", type = float, default = 0.001)
    # rate of sequencing error in simulation
    parser.add_argument("--sequencing_error_rate", type = float, default = 0.001)
    # whether to do junction rounding or not 
    parser.add_argument("--rounding", action = 'store_true')
    # The gender of the simulated sample
    parser.add_argument("--gender", default = 'Female')
    #
    main = argparse.ArgumentParser(prog = 'nebula', add_help = False)
    subparsers = main.add_subparsers(dest = 'command')
    ################## Preprocessing arguments
    gc_parser = subparsers.add_parser('misc', parents = [parser])
    unify_parser = subparsers.add_parser('unify', parents = [parser])
    genotype_parser = subparsers.add_parser('genotype', parents = [parser])
    simulation_parser = subparsers.add_parser('simulate', parents = [parser])
    preprocess_parser = subparsers.add_parser('preprocess', parents = [parser])
    ################## End of arguments
    args = main.parse_args()
    return args
