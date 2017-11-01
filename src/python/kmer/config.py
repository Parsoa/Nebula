import os
import pwd
import sys
import argparse

import colorama

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Configuration:

    kmer_cache_size = 10000

    class __impl:
        def __init__(self,
                        ksize,\
                        bed_file,\
                        genome_chm1,\
                        genome_hg19,\
                        genome_hg38,\
                        num_threads,\
                        khmer_num_tables,\
                        khmer_table_size,\
                        output_directory,\
                        sample_count_sevrer_port,\
                        sample_count_sevrer_port):
            self.ksize = ksize
            self.bed_file = bed_file
            self.genome_chm1 = genome_chm1
            self.genome_hg19 = genome_hg19
            self.genome_hg38 = genome_hg38
            self.num_threads = num_threads
            self.khmer_num_tables = khmer_num_tables
            self.khmer_table_size = khmer_table_size
            self.output_directory = output_directory
            self.sample_count_sevrer_port = sample_count_sevrer_port
            self.reference_count_sevrer_port = reference_count_sevrer_port

        def kmer_size(self):
            return self.ksize

    __instance = None

    def __init__(self, ksize=None, khmer_table_size=None, khmer_num_tables=None,\
            fastq_file=None, bed_file=None, output_directory=None, num_threads=None):
        if Configuration.__instance is None:
            Configuration.__instance = Configuration.__impl(ksize, khmer_table_size,\
                khmer_num_tables, fastq_file, bed_file, output_directory, num_threads)

    def __getattr__(self, attr):
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):
        return setattr(self.__instance, attr, value)

# ============================================================================================================================ #
# Configuration
# ============================================================================================================================ #

def configure():
    if sys.platform == "darwin":
        print('Running on Mac OS X')
        genome_chm1 = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1.samtoolsversion.head.small.fq'))
        genome_hg19 = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/hg19.fa'))
        genome_hg38 = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/hg38.fa'))
        num_threads = 1
        khmer_num_tables = 4
        khmer_table_size = 16e7
    else:
        print('Running on Linux')
        genome_chm1 = '/share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq'
        genome_hg19 = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref'
        genome_hg38 = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa'
        num_threads = 48
        khmer_num_tables = 4
        khmer_table_size = 16e9
    Configuration(
        ksize = 31,
        bed_file = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1_Lumpy.Del.100bp.bed')),
            # '../../../data/variations.bed')),
        genome_chm1 = genome_chm1,
        genome_hg19 = genome_hg19,
        genome_hg38 = genome_hg38,
        num_threads = num_threads,
        khmer_num_tables = khmer_num_tables,
        khmer_table_size = khmer_table_size,
        output_directory = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output')),
        sample_count_sevrer_port = 6985,
        reference_count_sevrer_port = 8569,
    )
    colorama.init()
 