from . import (
    bed,
    sets,
    config,
    commons,
    reference,
)

class StructuralVariation(object):

    def __init__(self, bed, radius):
        self.bed = bed
        self.radius = radius
        self.sequence = self.extract_base_sequence()

    def extract_base_sequence(self):
        track = copy.deepcopy(self.bed)
        track.start = track.start - self.radius
        track.end   = track.end   + self.radius
        sequence = bed.extract_sequence(track)

    def get_interval_boundaries(self, begin, end):