import copy

from . import (
    bed,
    sets,
    config,
    commons,
    reference,
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

    def get_interval_boundaries(self, begin, end):
        # test
        # track = copy.deepcopy(self.track)
        # track.start = track.start + begin
        # track.end   = track.end   + end
        # reference_boundaries, variation_boundaries = bed.extract_track_boundaries(track)
        #
        c = config.Configuration()
        # adjust boundaries
        # from [-R, R] to [0, 2R]
        begin = self.radius + begin
        # from [-R, R] to [2R, 0]
        end = self.radius - end 
        #
        l = len(self.sequence)
        seq = self.sequence[begin : l - end]
        # seq = seq[:c.ksize] + bed.complement_sequence((seq[c.ksize : -c.ksize])[::-1]) + seq[-c.ksize:]
        head = seq[0:2 * c.ksize]
        tail = seq[-2 * c.ksize:]
        # print(colorama.Fore.WHITE, variation_boundaries['head'])
        # print(colorama.Fore.BLUE, head)
        # print(colorama.Fore.WHITE, variation_boundaries['tail'])
        # print(colorama.Fore.BLUE, tail)
        # if variation_boundaries['head'] == head and variation_boundaries['tail'] == tail :
        #     print('match')
        return head, tail

