import copy
import time

from kmer import (
    bed,
    config
)

from kmer.kmers import *
from kmer.commons import *

import colorama

print('importing sv.py')

class SNP(object):

    def __init__(self, chrom, begin, end, variants):
        self.chrom = chrom
        self.begin = begin
        self.end = end
        self.variants = variants.strip().split('/')

class StructuralVariation(object):

    def __init__(self, track, radius):
        self.track = track
        self.radius = radius
        self.extract_base_sequence()

    def find_snps_within_boundaries(self, snps):
        print('here')
        pass

    # the end position itself is not included in the sequence
    def extract_base_sequence(self):
        c = config.Configuration()
        track = copy.deepcopy(self.track)
        # this is the largest sequence that we will ever need for this track
        # <- k bp -><- R bp -><-actual sequence-><- R bp -><- k bp ->
        track.start = track.start - self.radius - c.ksize
        track.end   = track.end   + self.radius + c.ksize
        self.sequence = bed.extract_sequence(track)
        #print(len(self.sequence), track.end - track.start + 1)

    def get_reference_signature_kmers(self, begin, end):
        c = config.Configuration()
        begin = (self.radius + c.ksize) + begin - c.ksize
        end = (len(self.sequence) - self.radius - c.ksize) + end + c.ksize
        seq = self.sequence[begin : end]
        #
        self.ref_head = seq[0:2 * c.ksize]
        self.ref_tail = seq[-2 * c.ksize:]
        kmers = extract_kmers(self.ref_head, self.ref_tail)
        return kmers

class Inversion(StructuralVariation):

    def get_signature_kmers(self, begin, end):
        c = config.Configuration()
        begin = (self.radius + c.ksize) + begin - c.ksize
        end = (len(self.sequence) - self.radius - c.ksize) + end + c.ksize
        seq = self.sequence[begin : end]
        # ends will overlap
        if begin >  end:
            return None, None
        #
        seq = seq[:c.ksize] + bed.complement_sequence((seq[c.ksize : -c.ksize])[::-1]) + seq[-c.ksize:]
        head = seq[0:2 * c.ksize]
        tail = seq[-2 * c.ksize:]
        # ends will overlap
        if 2 * c.ksize > len(seq) - 2 * c.ksize:
            return None, None
        #
        self.head = head
        self.tail = tail
        kmers = extract_kmers(head, tail)
        return kmers, {'head': head, 'tail': tail}

    def get_boundaries():
        return 

class Deletion(StructuralVariation):

    def find_snps_within_boundaries(self, snps):
        c = config.Configuration()
        begin = self.track.start - self.radius - c.ksize
        end = self.track.end + self.radius + c.ksize
        #pretty_print('checking for snps between', green(begin), 'and', green(end))
        events = {}
        for i in range(begin, end):
            if str(i) in snps and snps[str(i)].chrom == self.track.chrom:
                events[str(i)] = snps[str(i)]
        return events

    # <K bases><R bases>[self.track.start .....<self.track.end - 1>]<self.track.end><R bases><K bases>
    def get_signature_kmers_with_variation(self, event, begin, end):
        c = config.Configuration()
        begin = (self.radius + c.ksize) + begin - c.ksize
        end = (len(self.sequence) - self.radius - c.ksize) + end + c.ksize
        p = self.sequence[begin : end]
        p = p[:c.ksize] + p[-c.ksize:]
        # location in genome where the entire sequence for this sv begins
        offset = self.track.start - self.radius - c.ksize
        snp_index = int(event.begin) - offset
        # SNP falls inside the sequence for this track
        for variant in event.variants:
            if len(variant) == 1 and variant != '-':
                s = self.sequence[:snp_index] + variant + self.sequence[snp_index + 1:]
                if len(s) != len(self.sequence):
                    exit()
                s = s[begin : end]
                s = s[:c.ksize] + s[-c.ksize:]
                if s == p:
                    continue
                #print(self.track.chrom, ':', blue(self.track.start - self.radius - c.ksize) , '-', green(self.track.end + self.radius + c.ksize), ' @ ', cyan(event.begin), ',', variant, ' -> ', white(self.sequence[:snp_index]), cyan(variant), white(self.sequence[snp_index + 1:]), sep = '')
                kmers = extract_kmers(c.ksize, s)
                yield kmers, s, variant

    def get_signature_kmers(self, begin, end):
        #print('getting signature kmers')
        c = config.Configuration()
        begin = (self.radius + c.ksize) + begin - c.ksize
        end = (len(self.sequence) - self.radius - c.ksize) + end + c.ksize
        seq = self.sequence[begin : end]
        # ends will overlap
        if begin > end:
            return None, None
        #
        seq = seq[:c.ksize] + seq[-c.ksize:]
        kmers = extract_kmers(c.ksize, seq)
        return kmers, seq

