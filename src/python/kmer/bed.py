import io
import os
import pwd
import sys
import json
import time

from kmer import (
    config,
)

from kmer.commons import *

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

class BedTrack:

    def __init__(self, chrom, begin, end, **kwargs):
        self.chrom = chrom
        self.begin = begin
        self.end = end
        self.left_drift = 0
        self.right_drift = 0
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __str__(self):
        return self.chrom + '_' + str(self.begin) + '_' + str(self.end)

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

def track_from_name(name):
    tokens = name.lower().split('_')
    return BedTrack(tokens[0], int(tokens[1]), int(tokens[2]))

def load_tracks_from_file(path, keywords = []):
    tracks = []
    with open(path, 'r') as f:
        line = f.readline()
        while line:
            tokens = line.split()
            kwargs = {}
            for index, key in enumerate(keywords):
                if index + 3 < len(tokens):
                    kwargs[key] = tokens[3 + index]
            track = BedTrack(chrom = tokens[0], begin = int(tokens[1]), end = int(tokens[2]), **kwargs)
            tracks.append(track)
            line = f.readline()
    return tracks

