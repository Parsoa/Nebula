from __future__ import print_function

import copy
import time

from kmer import (
    bed,
    config
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import colorama

class SNP(object):

    def __init__(self, chrom, begin, end, variants):
        self.chrom = chrom
        self.begin = begin
        self.end = end
        self.variants = variants.strip().split('/')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class StructuralVariation(object):

    def __init__(self, track):
        self.track = track
        self.inner_kmers = None
        self.extract_base_sequence()

    def find_snps_within_boundaries(self, snps):
        pass

    # the begin position itself is not included in the sequence
    # the end position is included in the sequence
    # the origin will be -1, -1
    def extract_base_sequence(self):
        c = config.Configuration()
        track = copy.deepcopy(self.track)
        # this is the largest sequence that we will ever need for this track
        # <- k bp -><- R bp -><-actual sequence-><- R bp -><- k bp ->
        self.slack = c.insert_size - 2 * c.radius - c.ksize - 2 * c.read_length
        track.start = track.start - c.radius - c.ksize - c.read_length - self.slack
        track.end   = track.end   + c.radius + c.ksize + c.read_length + self.slack
        #self.sequence = bed.extract_sequence(track)
        #print(green(self.sequence))
        chromosome = extract_chromosome(track.chrom)
        self.sequence = chromosome[track.start - 1 : track.end - 1]
        #print(blue(self.sequence))

    # will return the same set of inner kmers for every breakpoint 
    def get_inner_kmers(self, counter, count, n):
        c = config.Configuration()
        begin = (c.radius + c.ksize + c.read_length + self.slack) + c.radius
        end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack) - c.radius
        if begin > end:
            return {}
        inner_seq = self.sequence[begin : end]
        # now count the kmers
        inner_kmers = extract_canonical_kmers(c.ksize, counter, count, inner_seq)
        items = inner_kmers.items()
        if len(inner_kmers) <= n:
            return inner_kmers
        else:
            return {kmer: inner_kmers[kmer] for kmer in list(map(lambda i: items[i][0], sorted(random.sample(range(0, len(inner_kmers)), n))))}

    def get_near_boundary_inner_kmers(self, counter = lambda x: 1, count = 1):
        c = config.Configuration()
        begin = (c.radius + c.ksize + c.read_length + self.slack)
        end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack)
        if begin > end:
            return {}
        inner_seq = self.sequence[begin : end]
        if end - begin < 6 * self.slack:
            return extract_kmers(c.ksize, counter, count, inner_seq)
        else:
            return extract_kmers(c.ksize, counter, count, inner_seq[: 3 * self.slack], inner_seq[-3 * self.slack :])

    # <L><Slack><K><R>|Event boundary|<R><Slack><L> ... <L><Slack><R>|Event Boundary|<R><K><Slack><L>
    def get_local_unique_kmers(self, counter):
        c = config.Configuration()
        left_end = self.sequence[:self.slack + c.read_length]
        right_end = self.sequence[-self.slack - c.read_length:]
        #print(green(self.sequence[:self.slack + c.read_length]) + cyan(self.sequence[self.slack + c.read_length: begin]) + white(self.sequence[begin : end]) + cyan(self.sequence[end : end + c.radius + c.ksize]) + green(self.sequence[end + c.radius + c.ksize :]))
        left_local_unique_kmers = extract_kmers(c.ksize, counter, 1, left_end)
        right_local_unique_kmers = extract_kmers(c.ksize, counter, 1, right_end)
        return right_local_unique_kmers, left_local_unique_kmers

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Inversion(StructuralVariation):

    def get_signature_kmers(self, begin, end, counter):
        c = config.Configuration()
        begin = (c.radius + c.ksize + c.read_length + self.slack) + begin - c.ksize
        end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack) + end + c.ksize
        seq = self.sequence[begin : end]
        # ends will overlap
        if begin >  end:
            return None, None
        #
        seq = seq[:c.ksize] + complement_sequence((seq[c.ksize : -c.ksize])[::-1]) + seq[-c.ksize:]
        head = seq[0:2 * c.ksize]
        tail = seq[-2 * c.ksize:]
        # ends will overlap
        if 2 * c.ksize > len(seq) - 2 * c.ksize:
            return None, None
        #
        self.head = head
        self.tail = tail
        kmers = extract_kmers(c.ksize, counter, 0, head, tail)
        return kmers, head + tail

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Deletion(StructuralVariation):

    def get_signature_kmers(self, begin, end, counter):
        #print('getting signature kmers')
        c = config.Configuration()
        begin = (c.radius + c.ksize + c.read_length + self.slack) + begin - c.ksize
        end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack) + end + c.ksize
        seq = self.sequence[begin : end]
        # ends will overlap
        if begin > end:
            return None, None
        #
        seq = seq[:c.ksize] + seq[-c.ksize:]
        boundary_kmers = extract_kmers(c.ksize, counter, 0, seq)
        return boundary_kmers, seq
