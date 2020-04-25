from __future__ import print_function

import io
import os
import pwd
import sys
import json
import time
import subprocess

from nebula import (
    config,
)

from nebula.kmers import *
from nebula.debug import *
from nebula.logger import *
from nebula.chromosomes import *

print = pretty_print

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

class BedTrack:

    def __init__(self, chrom, begin, end, fields = []):
        self.chrom = chrom
        self.begin = begin
        self.end = end
        self.fields = []
        for (k, v) in fields:
            self[k] = v
        if not hasattr(self, 'svtype'):
            self.svtype = 'DEL'
        if not hasattr(self, 'id'):
            self.id = self.svtype + '@' + self.chrom + '_' + str(self.begin) + '_' + str(self.end)

    def lift(self):
        with open('liftover_tmp.bed', 'w') as b:
            b.write(self.chrom + '\t' + str(self.begin) + '\t' + str(self.end) + '\n')
        command = '/home/pkhorsand/local/bin/liftOver ' + 'liftover_tmp.bed' + ' /afs/genomecenter.ucdavis.edu/home/pkhorsand/hg19ToHg38.over.chain liftover_res.bed liftover_un.bed'
        output = subprocess.call(command, shell = True)
        track = load_tracks_from_file('liftover_res.bed')[0]
        os.remove('liftover_tmp.bed')
        os.remove('liftover_res.bed')
        os.remove('liftover_un.bed')
        return track

    @staticmethod
    def json_serialize(self):
        d = {}
        d['end'] = self.end
        d['begin'] = self.begin
        d['chrom'] = self.chrom
        for (k, v) in self.fields:
            d[k] = v
        return d

    def serialize(self):
        s = ''
        s += str(self.chrom) + '\t' + str(self.begin) + '\t' + str(self.end)
        for key in self.fields:
            s += '\t' + str(self[key])
        s += '\n'
        return s

    def header(self):
        s = '#CHROM\tBEGIN\tEND'
        for key in self.fields:
            s += '\t' + key.upper()
        s += '\n'
        return s

    def __str__(self):
        return self.id

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        key = key.lower()
        setattr(self, key, value)
        if not key in self.fields:
            self.fields.append(key)

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

# sorts a dictionary of tracks into a list
def sort_tracks(tracks):
   _tracks = [tracks[track] for track in tracks]
   return sorted(sorted(_tracks, key = lambda x: x.begin), key = lambda y: y.chrom)

def filter_overlapping_tracks(tracks):
    i = 0
    remove = []
    tracks = sort_tracks(tracks)
    while i < len(tracks):
        for j in range(i + 1, len(tracks)):
            if tracks[j].chrom != tracks[i].chrom:
                i = j
                break
            if tracks[j].begin <= tracks[i].end:
                remove.append(j)
                user_print_warning(str(tracks[j]), 'overlaps', blue(str(tracks[i])))
                continue
            if tracks[j].begin - tracks[i].end < 1000:
                remove.append(j)
                user_print_warning(str(tracks[j]), 'is too close to', blue(str(tracks[i])))
                continue
            else:
                i = j
                break
        if j == len(tracks) - 1:
            break
    n = 0
    for index in sorted(remove):
        tracks.pop(index - n)
        n = n + 1
    return tracks

def track_from_id(name):
    svtype, coords = name.split('@')
    tokens = coords.split('_')
    return BedTrack(tokens[0], int(tokens[1]), int(tokens[2]), [('SVTYPE', svtype)])

# keywords is an array of tuples: (name, default, transformation)
def load_tracks_from_file(path, parse_header = True, keywords = []):
    tracks = []
    with open(path, 'r') as f:
        if parse_header:
            header = f.readline()
            fields = header.lower().split()
        line = f.readline()
        while line:
            tokens = line.split()
            kwargs = []
            if len(tokens) > 3:
                if parse_header:
                    for i in range(3, len(tokens)):
                        kwargs.append((fields[i], tokens[i]))
                else:
                    for index, pair in enumerate(keywords):
                        if index + 3 < len(tokens):
                            if len(pair) >= 3:
                                kwargs.append((pair[0], pair[2](tokens[3 + index])))
                            else:
                                kwargs.append((pair[0], tokens[3 + index]))
                        else:
                            kwargs.append((pair[0], pair[1]))
                track = BedTrack(chrom = tokens[0], begin = int(tokens[1]), end = int(tokens[2]), fields = kwargs)
                tracks.append(track)
            else:
                track = BedTrack(chrom = tokens[0], begin = int(tokens[1]), end = int(tokens[1]))
                tracks.append(track)
            line = f.readline()
    return tracks

def load_tracks_from_file_as_dict(path, parse_header = True, keywords = []):
    return {str(track): track for track in load_tracks_from_file(path, parse_header, keywords)}


