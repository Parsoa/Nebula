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
        self.sequence = chromosome[track.start: track.end]
        #print(blue(self.sequence))

    def get_inner_kmers(self, counter, count, n, begin = 0, end = 0, overlap = True, canonical = True):
        c = config.Configuration()
        begin = (c.radius + c.ksize + c.read_length + self.slack) + begin# + c.radius
        end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack) + end# - c.radius
        if begin > end:
            print(yellow('Event too short, no inner kmers exist'))
            return {}
        inner_seq = self.sequence[begin : end + c.ksize - 1]
        #print(inner_seq)
        #print(len(inner_seq))
        inner_kmers = c_extract_canonical_kmers(c.ksize, counter, count, overlap, inner_seq)
        #return inner_kmers
        if len(inner_kmers) <= n:
            return inner_kmers
        else:
            items = sorted(inner_kmers.items(), key = lambda item: item[1])[0:n]
            return {item[0]: item[1] for item in items}

    #def get_near_boundary_inner_kmers(self, counter = lambda x: 1, count = 1):
    #    c = config.Configuration()
    #    begin = (c.radius + c.ksize + c.read_length + self.slack)
    #    end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack)
    #    if begin > end:
    #        return {}
    #    inner_seq = self.sequence[begin : end]
    #    if end - begin < 6 * self.slack:
    #        return c_extract_kmers(c.ksize, counter, count, inner_seq)
    #    else:
    #        return c_extract_kmers(c.ksize, counter, count, inner_seq[: 3 * self.slack], inner_seq[-3 * self.slack :])

    ## <L><Slack><K><R>|Event boundary|<R><Slack><L> ... <L><Slack><R>|Event Boundary|<R><K><Slack><L>
    #def get_local_unique_kmers(self, counter):
    #    c = config.Configuration()
    #    left_end = self.sequence[:self.slack + c.read_length]
    #    right_end = self.sequence[-self.slack - c.read_length:]
    #    #print(green(self.sequence[:self.slack + c.read_length]) + cyan(self.sequence[self.slack + c.read_length: begin]) + white(self.sequence[begin : end]) + cyan(self.sequence[end : end + c.radius + c.ksize]) + green(self.sequence[end + c.radius + c.ksize :]))
    #    left_local_unique_kmers = extract_kmers(c.ksize, left_end)
    #    right_local_unique_kmers = extract_kmers(c.ksize, right_end)
    #    return right_local_unique_kmers, left_local_unique_kmers

    def get_boundary_kmers(self, begin, end, counter, count):
        pass

    def extract_boundary_gapped_kmers(self, counter = lambda x: 1, count = 1):
        c = config.Configuration()
        begin = (c.radius + c.ksize + c.read_length + self.slack)
        end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack)
        outer_gapped_kmers = {}
        inner_gapped_kmers = {}
        g = c.gap / 2
        k = (c.ksize / 2) * 2
        b = self.sequence[begin - c.ksize: begin + c.ksize]
        for i in range(0, c.gap):
            kmer = self.sequence[begin - i - 15: begin - i + 5 + 15]
            #kmer = kmer[:15] + kmer[20:]
            inner_gapped_kmers[kmer] = True
        e = self.sequence[end - c.ksize: end + c.ksize]
        for i in range(0, c.gap):
            kmer = self.sequence[end - 1 - i - 15: end - 1 - i + 5 + 15]
            #kmer = kmer[:15] + kmer[20:]
            inner_gapped_kmers[kmer] = True
        seq = self.sequence[: begin] + self.sequence[end :]
        base = c.radius + c.ksize + c.read_length + self.slack
        for i in range(0, c.gap):
            kmer = seq[base - i - 15: base - i + 5 + 15]
            #kmer = kmer[:15] + kmer[20:]
            outer_gapped_kmers[kmer] = True
        return {'inner': inner_gapped_kmers, 'outer': outer_gapped_kmers, 'begin': b, 'end': e}

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
        kmers = extract_kmers(counter, 0, head, tail)
        return kmers, head + tail

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Deletion(StructuralVariation):

    def get_boundary_kmers(self, begin, end, counter, count = 0):
        #print('getting signature kmers')
        c = config.Configuration()
        begin = (c.radius + c.ksize + c.read_length + self.slack) + begin - c.ksize
        end = (len(self.sequence) - c.radius - c.ksize - c.read_length - self.slack) + end + c.ksize
        seq = self.sequence[begin : end].upper()
        # ends will overlap
        if begin > end:
            return None, None
        #
        #print(green(seq[:c.ksize]), white(seq[c.ksize:-c.ksize]), green(seq[-c.ksize:]))
        seq = seq[:c.ksize] + seq[-c.ksize:]
        boundary_kmers = c_extract_canonical_kmers(c.ksize, counter, 0, True, seq)
        return boundary_kmers, seq
