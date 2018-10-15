from __future__ import print_function

import copy
import time

from kmer import (
    config
)

from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

import colorama

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class StructuralVariation(object):

    def __init__(self, track):
        self.chrom = track.chrom
        self.begin = track.begin
        self.end = track.end
        self.inner_kmers = None
        self.extract_base_sequence()

    # the begin position itself is not included in the sequence
    # the end position is included in the sequence
    # the origin will be -1, -1
    def extract_base_sequence(self):
        c = config.Configuration()
        # this is the largest sequence that we will ever need for this track
        # <- k bp -><- R bp -><-actual sequence-><- R bp -><- k bp ->
        self.slack = c.insert_size - 2 * c.radius - c.ksize - 2 * c.read_length
        begin = self.begin - c.radius - c.ksize - c.read_length - self.slack
        end = self.end   + c.radius + c.ksize + c.read_length + self.slack
        chromosome = extract_chromosome(self.chrom)
        self.sequence = chromosome[begin: end]

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
        print('getting inner kmers')
        inner_kmers = c_extract_canonical_kmers(c.ksize, counter, count, overlap, inner_seq)
        #return inner_kmers
        if len(inner_kmers) <= n:
            return inner_kmers
        else:
            items = sorted(inner_kmers.items(), key = lambda item: item[1])[0:n]
            return {item[0]: item[1] for item in items}

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
        # each half is 15 bases plus 10 bases in between, say nearly 45: 25 bases remain for each end
        b = self.sequence[begin - 15 - 2 - 25: begin + 3 + 15 + 25]
        kmer = self.sequence[begin - 15 - 2: begin + 3 + 15]
        prefix = self.sequence[begin - 15 - 2 - 25: begin - 15 - 2]
        suffix = self.sequence[begin + 3 + 15: begin + 3 + 15 + 25]
        inner_gapped_kmers[kmer] = {'indicators': self.generate_kmer_mask(prefix, suffix), 'prefix': prefix, 'suffix': suffix}
        #
        e = self.sequence[end - 15 - 2 - 25: end + 3 + 15 + 25]
        kmer = self.sequence[end - 15 - 2: end + 3 + 15]
        prefix = self.sequence[end - 15 - 2 - 25: end - 15 - 2]
        suffix = self.sequence[end + 3 + 15: end + 3 + 15 + 25]
        inner_gapped_kmers[kmer] = {'indicators': self.generate_kmer_mask(prefix, suffix), 'prefix': prefix, 'suffix': suffix}
        #
        kmer = self.sequence[begin - 2 - 15: begin + 3] + self.sequence[end - 2: end + 3 + 15]
        prefix = self.sequence[begin - 15 - 2 - 25: begin - 15 - 2]
        suffix = self.sequence[end + 3 + 15: end + 3 + 15 + 25]
        outer_gapped_kmers[kmer] = {'indicators': self.generate_kmer_mask(prefix, suffix), 'prefix': prefix, 'suffix': suffix}
        return {'inner': inner_gapped_kmers, 'outer': outer_gapped_kmers, 'begin': b, 'end': e, 'sequence': self.sequence}

    def generate_kmer_mask(self, left, right):
        c = config.Configuration()
        masks = {}
        for seq in [left, right]:
            for j in range(0, 5):
                indices = []
                while len(indices) != 21:
                    i = random.randint(0, len(seq) - 1)
                    if not i in indices:
                        indices.append(i)
                indices = sorted(indices)
                indices.insert(0, '')
                mask = reduce(lambda x, y: x + seq[y], indices)
                masks[mask] = True
        return masks

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
