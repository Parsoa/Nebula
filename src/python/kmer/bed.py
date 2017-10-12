import os
import pwd
import sys
import time

from . import (
    sets,
    config,
    reference
)

import khmer
import pybedtools

def measure_time(f):
    def wrapper():
        start = time.clock()
        result = f()
        end = time.clock()
        print('took ', end - start)
        return result
    return wrapper

def conditionally(f):
    print(__name__)
    if __name__ == '__main__':
        print('not profiling')
        return f
    print('profiling')
    return profile(f)

@conditionally
def read_tracks_from_bed_file():
    c = config.Configuration()
    bedtools = pybedtools.BedTool(c.bed_file)
    # 
    # reference_counttable = count_reference_kmers()
    # sample_counttable = count_sample_kmers()
    #
    for track in bedtools:
        print('--------------------------------------------------------')
        print('track: ', str(track).strip())
        # 
        reference_boundaries, variation_boundaries = extract_track_boundaries(track)
        reference_counttable, reference_kmers = extract_kmers_of_interest(reference_boundaries)
        variation_counttable, variation_kmers = extract_kmers_of_interest(variation_boundaries)
        #
        print('reference segment: ', len(reference_kmers.keys()), ' @ ',\
            reference_boundaries['head'], ' ... ', reference_boundaries['tail'])
        print('variation segment: ', len(variation_kmers.keys()), ' @ ',\
            variation_boundaries['head'], ' ... ', variation_boundaries['tail'])
        # 
        intersection = sets.calc_dictionary_intersection(reference_kmers, variation_kmers)
        print('reference/variation intersection: ', len(intersection), ' kmers')
        if intersection:
            print(intersection)
        #
        # reference_score = len(calc_similarity_score(
        #     sets.calc_dictionary_difference(reference_kmers, variation_kmers), sample_counttable))
        # print('fastq/reference similarity: ', reference_score)
        # variation_score = len(calc_similarity_score(
        #     sets.calc_dictionary_difference(variation_kmers, reference_kmers), sample_counttable))
        # print('fastq/variation similarity: ', variation_score)
        # print('========================================')
        # print('decision: ', ('reference' if reference_score > variation_score else
        #             ('variation' if reference_score < variation_score else 'undecisive')))

def extract_kmers_of_interest(boundaries):
    c = config.Configuration()
    # 
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    #
    kmers = {}
    #
    counttable.consume(boundaries['head'])
    counttable.consume(boundaries['tail'])
    #
    for sequence in [boundaries['head'], boundaries['tail']]:
        for kmer in counttable.get_kmers(sequence) :
            kmers[kmer] = counttable.get_kmer_counts(kmer)
    return counttable, kmers

def count_kmers_from_file(file):
    c = config.Configuration()
    #
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counttable.consume_seqfile(file)
    #
    return counttable, nkmers

@measure_time
def count_reference_kmers():
    print('reading reference genome ... ')
    c = config.Configuration()
    #
    counttable, nkmers = count_kmers_from_file(reference.ReferenceGenome().path)
    #
    print('reference counttable cached: ', nkmers, ' kmers recorded')
    return counttable

@measure_time
def count_sample_kmers():
    print('reading sample genome ...')
    c = config.Configuration()
    # 
    print('looking for cached counttable')
    cache = c.fastq_file + '.ct'
    if os.path.isfile(cache):
        print('found at ', cache)
        return khmer.Counttable.load(cache)
    #
    print('not found, creating counttable ' + cache)
    counttable, nkmers = count_kmers_from_file(c.fastq_file)
    counttable.save(cache)
    #
    print('sample counttable cached: ', nkmers, ' kmers recorded')
    return counttable

def count_variation_kmers(reference_counttable, track):
    #TODO implement
    sequence = extract_reference_sequence(track)

def extract_track_boundaries(track):
    c = config.Configuration()
    interval = pybedtools.Interval(chrom = track.chrom, start = track.start - c.ksize,\
        end = track.end + c.ksize)
    bedtool = pybedtools.BedTool(str(interval), from_string = True)
    f = open((bedtool.sequence(fi = reference.ReferenceGenome().fasta)).seqfn)
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

def calc_similarity_score(kmers, counttable):
    print('========================================')
    result = {}
    for kmer in kmers:
        if counttable.get_kmer_counts(kmer)[0] != 0 :
            print(kmer, 'sample: ', '{:04d}'.format(counttable.get_kmer_counts(kmer)[0]))
            result[kmer] = True
    return result

def configure():
    if sys.platform == "darwin":
        print('Running on Mac OS X')
        khmer_table_size = 16e7
        khmer_num_tables = 4
        reference.ReferenceGenome(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/hg38.fa')))
        fastq_file = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1.samtoolsversion.head.tiny.fq'))
    else:
        print('Running on Linux')
        khmer_table_size = 16e9
        khmer_num_tables = 4
        fastq_file = '/share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq'
        reference.ReferenceGenome('/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa')
    config.Configuration(
        ksize = 25,
        khmer_table_size = khmer_table_size,
        khmer_num_tables = khmer_num_tables,
        fastq_file = fastq_file,
        bed_file = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1.inversions_hg38.bed'))
    )

if __name__ == '__main__':
    configure()
    read_tracks_from_bed_file()
 