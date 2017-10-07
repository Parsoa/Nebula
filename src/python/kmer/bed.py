import sys
import os
import time
import pwd

from . import (
    reference,
    config,
    sets,
)

import pybedtools
import khmer

# @profile
def read_tracks_from_bed_file(path):
    c = config.Configuration()
    bedtools = pybedtools.BedTool(path)
    # 
    # reference_counttable = count_reference_kmers()
    # sample_counttable = count_sample_kmers()
    #
    for track in bedtools:
        print('track: ', track)
        reference_boundaries, variation_boundaries = extract_track_boundaries(track)
        # 
        reference_counttable, reference_kmers = extract_kmers_of_interest(reference_boundaries)
        variation_counttable, variation_kmers = extract_kmers_of_interest(variation_boundaries)
        #
        print('reference kmers: ', len(reference_kmers.keys()))
        print('variation kmers: ', len(variation_kmers.keys()))
        print('reference/variation intersection: ', len(sets.calc_dictionary_intersection(reference_kmers, variation_kmers)))
        #
        reference_score = len(calc_jaccard_similarity(
            sets.calc_dictionary_difference(kmers, inverse_kmers), fastq_counts))
        print('fastq/reference similarity: ', reference_score)
        inverse_score = len(calc_jaccard_similarity(
            sets.calc_dictionary_difference(inverse_kmers, kmers), fastq_counts))
        print('fastq/inverse similarity: ', inverse_score)
        print('decision: ', ('reference' if reference_score > inverse_score else ('inverse' if reference_score < inverse_score else 'undecisive')))

def extract_kmers_of_interest(boundaries):
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    #
    kmers = {}
    #
    countt.consume(boundaries['head'])
    countt.consume(boundaries['tail'])
    #
    for sequence in [boundaries['head'], boundaries['tail']]:
        for kmer in counttable.get_kmers(sequence) :
            kmers[kmer] = countt.get_kmer_counts(kmer)
    return counttable, kmers

def count_kmers_from_file(file):
    c = config.Configuration()
    #
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counttable.consume_seqfile(file)
    #
    return counttable, nkmer

def count_reference_kmers():
    print('reading reference genome ... ')
    c = config.Configuration()
    start = time.clock()
    #
    counttable, nkmers = count_kmers_from_file(reference.ReferenceGenome().fasta)
    #
    end = time.clock()
    print('took ', end - start)
    print('kmers in reference genome: ', nkmers)
    return counttable, nkmers

def count_sample_kmers():
    print('reading sample ...')
    c = config.Configuration()
    start = time.clock()
    #
    counttable, nkmers = count_kmers_from_file(c.fastq_file)
    #
    end = time.clock()
    print('took ', end - start)
    print('sample cached, ', nkmers, ' kmers')

def count_variation_kmers(reference_counttable, track):
    #TODO implement
    sequence = extract_reference_sequence(track)

def extract_track_boundaries(track):
    c = config.Configuration()
    interval = pybedtools.Interval(chrom=track.chrom, start=track.start - c.ksize, end=track.end + c.ksize)
    bedtool = pybedtools.BedTool(str(interval), from_string=True)
    print('reference genome: ' + reference.ReferenceGenome().path)
    f = open(bedtool.sequence(fi = reference.ReferenceGenome().fasta))
    for i, line in enumerate(f):
        line = line.strip()
        if i == 1:
            # print('#', line, '#')
            head = line[0:2 * c.ksize]
            tail = line[-2 * c.ksize:]
            # print('head: ', head)
            # print('tail: ', tail)
            # inverse this sequence
            line = line[::-1]
            # print('#', line, '#')
            inverse_head = line[0:2 * c.ksize]
            inverse_tail = line[-2 * c.ksize:]
            # print('inverse head: ', inverse_head)
            # print('inverse tail: ', inverse_tail)
            return {'head': head.upper(), 'tail': tail.upper()}, {'head': inverse_head.upper(), 'tail': inverse_tail.upper()}

def calc_jaccard_similarity(kmers, countgraph):
    result = {}
    for kmer in kmers:
        if countgraph.get_kmer_counts(kmer)[0] != 0 :
            result[kmer] = True
    return result


def configure():
    # print('cwd: {}'.format(os.getcwd()))
    khmer_table_size = 5e8
    khmer_num_tables = 4
    if sys.platform == "darwin":
        print('Running on Mac OS X')
        reference.ReferenceGenome(os.path.join(os.path.dirname(__file__), '../../../data/hg38.fa'))
        fastq_file = os.path.join(os.path.dirname(__file__), '../../../data/CHM1.samtoolsversion.head.fq')
    else:
        print('Running on Linux')
        fastq_file = '/share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.head.fq'
        reference.ReferenceGenome('/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa')
    config.Configuration(
        ksize = 25,
        khmer_table_size = khmer_table_size,
        khmer_num_tables = khmer_num_tables,
        fastq_file = fastq_file
    )
    bed_file = os.path.join(os.path.dirname(__file__),'../../../data/CHM1.inversions_hg38.bed')
    read_tracks_from_bed_file(path = bed_file)

if __name__ == '__main__':
    configure()
 