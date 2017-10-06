import sys
import time

from kmer import (
    reference,
    config,
    sets,
)

import pybedtools
import khmer

def read_tracks_from_bed_file(path):
    c = config.Configuration()
    bedtools = pybedtools.BedTool(path)
    #
    print('reading sample ...')
    start = time.clock()
    fastq_counts = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = fastq_counts.consume_seqfile(c.fastq_file)
    end = time.clock()
    print('took ', end - start)
    print('sample cached, ', nkmers, ' kmers')
    #
    for track in bedtools:
        print('track: ', track)
        counts = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
        inverse_counts = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
        #
        kmers = {}
        inverse_kmers = {}
        #
        sequence = extract_reference_sequence(track)
        (head, tail), (inverse_head, inverse_tail) = extract_sequence_boundaries(sequence)
        print('reference: [', head, '.....', tail, ']')
        print('inverse: [', inverse_head, '.....', inverse_tail, ']')
        #
        counts.consume(head)
        counts.consume(tail)
        #
        inverse_counts.consume(inverse_head)
        inverse_counts.consume(inverse_tail)
        #
        for seq in [head, tail] :
            for kmer in counts.get_kmers(seq) :
                # print(kmer, ':', counts.get_kmer_counts(kmer))
                kmers[kmer] = True
        #
        for seq in [inverse_head, inverse_tail] :
            for kmer in inverse_counts.get_kmers(seq) :
                inverse_kmers[kmer] = True
        #
        print('reference kmers: ', len(kmers.keys()))
        print('inverse kmers: ', len(inverse_kmers.keys()))
        print('inverse/reference intersection: ', len(sets.calc_dictionary_intersection(kmers, inverse_kmers)))
        #
        reference_score = len(calc_jaccard_similarity(
            sets.calc_dictionary_difference(kmers, inverse_kmers), fastq_counts))
        print('fastq/reference similarity: ', reference_score)
        inverse_score = len(calc_jaccard_similarity(
            sets.calc_dictionary_difference(inverse_kmers, kmers), fastq_counts))
        print('fastq/inverse similarity: ', inverse_score)
        print('decision: ', ('reference' if reference_score > inverse_score else ('inverse' if reference_score < inverse_score else 'undecisive')))

    # print('common:', sets.print_dictionary_keys(sets.calc_dictionary_intersection(kmers, inverse_kmers)))

def extract_reference_sequence(track):
    # TODO: this might actually cause problem if it surpasses the boundaris of the chromosome
    c = config.Configuration()
    interval = pybedtools.Interval(chrom=track.chrom, start=track.start - c.ksize, end=track.end + c.ksize)
    # print('Interval: ', interval)
    bedtool = pybedtools.BedTool(str(interval), from_string=True)
    sequence = bedtool.sequence(fi=reference.ReferenceGenome().fasta)
    return sequence.seqfn

def extract_sequence_boundaries(sequence):
    # print(sequence)
    c = config.Configuration()
    f = open(sequence)
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
            return (head.upper(), tail.upper()), (inverse_head.upper(), inverse_tail.upper())

def calc_jaccard_similarity(kmers, countgraph):
    result = {}
    for kmer in kmers:
        if countgraph.get_kmer_counts(kmer)[0] != 0 :
            result[kmer] = True
    return result

def extract_sample_sequence(track):
    bedtool = pybedtools.BedTool(str(track), from_string=True)
    sequence = bedtool.sequence(fi=reference.ReferenceGenome().fasta)
    return sequence.seqfn

def count_kmers(sequence):
    c = config.Configuration()
    counts = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counts.consume_seqfile(sequence)
    print(nseqs, nkmers)
