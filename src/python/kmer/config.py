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
                        max_threads,\
                        variation_type,\
                        reference_genome,\
                        khmer_num_tables,\
                        khmer_table_size,\
                        output_directory,\
                        sample_count_server_port,\
                        reference_count_server_port):
            self.ksize = ksize
            self.bed_file = bed_file
            self.genome_chm1 = genome_chm1
            self.genome_hg19 = genome_hg19
            self.genome_hg38 = genome_hg38
            self.max_threads = max_threads
            self.variation_type = variation_type
            self.reference_genome = reference_genome
            self.khmer_num_tables = khmer_num_tables
            self.khmer_table_size = khmer_table_size
            self.output_directory = output_directory
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

def configure(reference_genome = 'hg38', bed_file = "CHM1_Lumpy.Del.100bp.bed", variation_type = 'Deletion'):
    if sys.platform == "darwin":
        print('Running on Mac OS X')
        genome_chm1 = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1.samtoolsversion.head.small.fq'))
        genome_hg19 = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/hg38.fa'))
        genome_hg38 = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/hg38.fa'))
        max_threads = 1
        khmer_num_tables = 4
        khmer_table_size = 16e7
    else:
        print('Running on Linux')
        genome_chm1 = '/share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq'
        genome_hg19 = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref'
        genome_hg38 = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa'
        max_threads = 48
        khmer_num_tables = 4
        khmer_table_size = 16e9
    #
    Configuration(
        ksize = 31,
        bed_file = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/' + bed_file)),
        genome_chm1 = genome_chm1,
        genome_hg19 = genome_hg19,
        genome_hg38 = genome_hg38,
        max_threads = max_threads,
        variation_type = variation_type,
        reference_genome = genome_hg19 if reference_genome == 'hg19' else genome_hg38,
        khmer_num_tables = khmer_num_tables,
        khmer_table_size = khmer_table_size,
        output_directory = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output')),
        sample_count_server_port = 6985,
        reference_count_server_port = 8569
    )
    colorama.init()
 