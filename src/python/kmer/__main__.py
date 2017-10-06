import os
import pwd
from sys import platform

from kmer import (
    bed,
    reference,
    config
)

if __name__ == '__main__':
    print('cwd: {}'.format(os.getcwd()))
    khmer_table_size = 5e8
    khmer_num_tables = 4
    if platform == "darwin":
        print('Running on Mac OS X')
        reference_genome = '/Users/' + pwd.getpwuid(os.getuid()).pw_name + '/Desktop/Davis/Projects/NebulousSerendipity/data/hg38.fa'
        fastq_file = '/Users/Parsoa/Desktop/Davis/Projects/NebulousSerendipity/data/CHM1.samtoolsversion.head.tiny.fq'
    else:
        print('Running on Linux')
        reference_genome = '/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa'
        fastq_file='/share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.head.tiny.fq'
    config.Configuration(
        ksize=25,
        khmer_table_size=khmer_table_size,
        khmer_num_tables=khmer_num_tables,
        fastq_file=fastq_file
    )
    bed_file = os.path.abspath('../../data/CHM1.inversions_hg38.bed')
    reference.ReferenceGenome(path=reference_genome)
    bed.read_tracks_from_bed_file(path=bed_file)