import copy

from kmer import (
    bed,
    sets,
    config,
    commons,
)

import colorama

class StructuralVariation(object):

    def __init__(self, track, radius):
        self.track = track
        self.radius = radius
        self.extract_base_sequence()

    def extract_base_sequence(self):
        c = config.Configuration()
        track = copy.deepcopy(self.track)
        # this is the largest sequence that we will ever need for this track
        # <- k bp -><- R bp -><-actual sequence-><- R bp -><- k bp ->
        track.start = track.start - self.radius - c.ksize
        track.end   = track.end   + self.radius + c.ksize
        self.sequence = bed.extract_sequence(track)

    def get_reference_signature_kmers(self, begin, end):
        c = config.Configuration()
        # adjust from [-R, R] to [0, 2R]
        begin = self.radius + begin
        # adjust from [-R, R] to [2R, 0]
        end = self.radius - end 
        #
        seq = self.sequence[begin : len(self.sequence) - end]
        self.ref_head = seq[0:2 * c.ksize]
        self.ref_tail = seq[-2 * c.ksize:]
        kmers = count_server.count_kmers_exact_list(self.ref_head, self.ref_tail)
        return kmers

class Inversion(StructuralVariation):

    def get_signature_kmers(self, begin, end, complement):
        c = config.Configuration()
        # adjust from [-R, R] to [0, 2R]
        begin = self.radius + begin
        # adjust from [-R, R] to [2R, 0]
        end = self.radius - end 
        #
        seq = self.sequence[begin : len(self.sequence) - end]
        if complement:
            seq = seq[:c.ksize] + bed.complement_sequence((seq[c.ksize : -c.ksize])[::-1]) + seq[-c.ksize:]
        head = seq[0:2 * c.ksize]
        tail = seq[-2 * c.ksize:]
        # ends will overlap
        if 2 * c.ksize > len(seq) - 2 * c.ksize:
            return None, None
        self.head = head
        self.tail = tail
        kmers = count_server.count_kmers_exact_list(head, tail)
        return kmers, {'head': head, 'tail': tail}

    def get_boundaries():
        return 

class Deletion(StructuralVariation):

    def get_signature_kmers(self, begin, end, delete):
        c = config.Configuration()
        # adjust from [-R, R] to [0, 2R]
        begin = self.radius + begin
        # adjust from [-R, R] to [2R, 0]
        end = self.radius - end 
        #
        seq = self.sequence[begin : len(self.sequence) - end]
        if delete:
            seq = seq[:c.ksize] + seq[len(self.sequence) - c.ksize:]
        kmers = count_server.count_kmers_exact_list(seq)
        return kmers, seq

