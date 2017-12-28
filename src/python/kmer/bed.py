import io
import os
import pwd
import sys
import json
import time

from . import (
    sets,
    config,
    commons,
)

import khmer
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

def export_tracks():
    c = config.Configuration()
    bedtools = pybedtools.BedTool(c.bed_file)
    #
    kmers = {}
    for track in bedtools:
        print(colorama.Fore.GREEN + '--------------------------------------------------------')
        #
        key = track.chrom.strip() + '_' + str(track.start).strip()  + '_'+ str(track.end).strip()
        print(colorama.Fore.GREEN + 'track: ', key)
        # 
        reference_boundaries, variation_boundaries = extract_track_boundaries(track)
        reference_kmers = count_boundary_kmers(reference_boundaries)
        variation_kmers = count_boundary_kmers(variation_boundaries)
        #
        kmers[key] = {
            'start': track.start,
            'end'  : track.end,
            'reference': {
                'head':  reference_boundaries['head'],
                'tail':  reference_boundaries['tail'],
                'kmers': reference_kmers
            },
            'variation': {
                'head':  variation_boundaries['head'],
                'tail':  variation_boundaries['tail'],
                'kmers': variation_kmers
            }
        }
    # print(colorama.Fore.GREEN, kmers)
    path = os.path.join(c.output_directory, 'boundaries.json')
    print(colorama.Fore.GREEN + 'exporting tracks to ', path)
    with open(path, 'w') as json_file:
        json.dump(kmers, json_file, indent = 4, separators = (',', ': '))

def import_tracks():
    c = config.Configuration()
    with io.open(os.path.join(c.output_directory, 'boundaries.json'), 'r') as json_file:
        tracks = json.load(json_file)
        return tracks

def extract_sequence(track):
    c = config.Configuration()
    interval = pybedtools.Interval(chrom = track.chrom, start = track.start,\
        end = track.end)
    bedtool = pybedtools.BedTool(str(interval), from_string = True)
    f = open(bedtool.sequence(c.reference_genome).seqfn)
    for i, line in enumerate(f):
        if i == 1:
            line = line.strip().upper()
            return line

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
