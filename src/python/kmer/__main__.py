import os
import pwd

from kmer import (
    bed,
    reference,
    config,
)

if __name__ == '__main__':
    khmer_table_size = 5e8
    khmer_num_tables = 4
    config.Configuration(
        ksize=25,
        khmer_table_size=khmer_table_size,
        khmer_num_tables=khmer_num_tables,
        fastq_file='/share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq'
    )
    reference.ReferenceGenome(path='/Users/' + pwd.getpwuid(os.getuid()).pw_name + '/Desktop/Davis/Projects/NebulousSerendipity/data/hg38.fa')
    bed.read_tracks_from_bed_file(path='/Users/' + pwd.getpwuid(os.getuid()).pw_name + '/Desktop/Davis/Projects/NebulousSerendipity/data/CHM1.inversions_hg38.bed')
