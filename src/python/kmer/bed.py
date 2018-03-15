import io
import os
import pwd
import sys
import json
import time

from kmer import (
    sets,
    config,
    commons,
)

import colorama
import pybedtools

print('importing bed.py')
# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

class BedTrack:

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = chrom + ':' + str(start) + '-' + str(end)

    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.__dict__ == other.__dict__
        return False

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

def read_tracks(path):
    c = config.Configuration()
    bedtools = pybedtools.BedTool(path)
    tracks = {}
    for track in bedtools:
        bed_track = BedTrack(chrom = track.chrom, start = track.start, end = track.end)
        if bed_track not in tracks:
            tracks[bed_track] = True
    return tracks

def extract_sequence(track):
    c = config.Configuration()
    interval = pybedtools.Interval(chrom = track.chrom, start = track.start,\
        end = track.end)
    bedtool = pybedtools.BedTool(str(interval), from_string = True)
    f = open(bedtool.sequence(c.reference_genome).seqfn)
    sequence = ''
    for i, line in enumerate(f):
        if i >= 1:
            line = line.strip().upper()
            sequence += line
    return sequence

def parse_bed_file(f):
    line = f.readline()
    while line:
        tokens = line.split()
        yield tokens
        line = f.readline()

def reverse_complement_sequence(seq):
    return complement_sequence(seq[::-1])

def complement_sequence(seq):
    # A-> C and C->A
    seq = seq.replace('A', 'Z')
    seq = seq.replace('T', 'A')
    seq = seq.replace('Z', 'T')
    #
    seq = seq.replace('G', 'Z')
    seq = seq.replace('C', 'G')
    seq = seq.replace('Z', 'C')
    #
    return seq
