import sys

from kmer import (
    reference,
    config,
)

import pybedtools
import khmer

def read_tracks_from_bed_file(path):
    c = config.Configuration()
    counts = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    bedtools = pybedtools.BedTool(path)
    for track in bedtools:
        sequence = extract_reference_sequence(track)
        head, tail = extract_sequence_boundaries(sequence)
        print(counts.get_kmers(head))
        counts.consume(head)
        counts.consume(tail)

def extract_reference_sequence(track):
    # TODO: this might actually cause problem if it surpasses the boundaris of the chromosome
    c = config.Configuration()
    interval = pybedtools.Interval(chrom=track.chrom, start=track.start - c.ksize, end=track.end + c.ksize)
    # print('Interval: ', interval)
    bedtool = pybedtools.BedTool(str(interval), from_string=True)
    sequence = bedtool.sequence(fi=reference.ReferenceGenome().fasta)
    return sequence.seqfn

def extract_sequence_boundaries(sequence):
    print(sequence)
    c = config.Configuration()
    f = open(sequence)
    for i, line in enumerate(f):
        line = line.strip()
        if i == 1:
            head = line[0:2 * c.ksize]
            tail = line[-2 * c.ksize:]
            print('head: ', head)
            print('tail: ', tail)
            return head, tail

def extract_sample_sequence(track):
    bedtool = pybedtools.BedTool(str(track), from_string=True)
    sequence = bedtool.sequence(fi=reference.ReferenceGenome().fasta)
    # print(open(sequence.seqfn).read())
    return sequence.seqfn

def count_kmers(sequence):
    c = config.Configuration()
    counts = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counts.consume_seqfile(sequence)
    print(nseqs, nkmers)
