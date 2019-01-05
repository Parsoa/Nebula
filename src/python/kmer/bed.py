import io
import os
import pwd
import sys
import json
import time

from kmer import (
    config,
)

from kmer.debug import *
from kmer.commons import *

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

class BedTrack:

    def __init__(self, chrom, begin, end, **kwargs):
        self.chrom = chrom
        self.begin = begin
        self.end = end
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.kwargs = kwargs

    def __str__(self):
        return self.chrom + '_' + str(self.begin) + '_' + str(self.end)

    def export(self):
        s = self.chrom + '\t' + str(self.begin) + '\t' + str(self.end)
        for k, v in enumerate(self.kwargs):
            s += '\t' + str(self.kwargs[v])
        s += '\n'
        return s

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

def track_from_name(name):
    tokens = name.split('_')
    return BedTrack(tokens[0], int(tokens[1]), int(tokens[2]))

# keywords is an array of tuples: (name, default, transformation)
def load_tracks_from_file(path, keywords = []):
    tracks = []
    with open(path, 'r') as f:
        line = f.readline()
        while line:
            tokens = line.split()
            kwargs = {}
            if len(tokens) >= 3:
                for index, pair in enumerate(keywords):
                    if index + 3 < len(tokens):
                        if len(pair) >= 3:
                            kwargs[pair[0]] = pair[2](tokens[3 + index])
                        else:
                            kwargs[pair[0]] = tokens[3 + index]
                    else:
                        kwargs[pair[0]] = pair[1]
                track = BedTrack(chrom = tokens[0], begin = int(tokens[1]), end = int(tokens[2]), **kwargs)
                tracks.append(track)
            else:
                track = BedTrack(chrom = tokens[0], begin = int(tokens[1]), end = int(tokens[1]))
                tracks.append(track)
            line = f.readline()
    return tracks

def load_tracks_from_file_as_dict(path, keywords = []):
    return {str(track): track for track in load_tracks_from_file(path, keywords)}
