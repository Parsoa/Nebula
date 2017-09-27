import sys

from kmer import (
    reference,
    config,
)

import pybedtools
import khmer

def read_tracks_from_bed_file(path):
    bedtools = pybedtools.BedTool(path)
    for track in bedtools:
        sequence = extract_reference_sequence(track)
        count_kmers(sequence)

def extract_reference_sequence(track):
    bedtool = pybedtools.BedTool(str(track), from_gok5443w;'rrfnstring=True)
    sequence = bedtool.sequence(fi=reference.ReferenceGenome().fasta)
    # print(open(sequence.seqfn).read())
    return sequence.seqfn

def count_kmers(sequence):
    c = config.Configuration()
    counts = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counts.consume_seqfile(sequence)
    print(nseqs, nkmers)
