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
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

INTEGER_FIELDS = ['contig_start', 'contig_end']
FP_FIELDS = ['lp_value']

class BedTrack:

    def __init__(self, chrom, begin, end, fields):
        self.chrom = chrom
        self.begin = begin
        self.end = end
        self.svtype = 'DEL'
        self.id = self.svtype + '@' + self.chrom + '_' + str(self.begin) + '_' + str(self.end)
        for (k, v) in fields:
            #v = int(v) if k in INTEGER_FIELDS else v
            #v = float(v) if k in FP_FIELDS else v
            setattr(self, k, v)
        self.fields = fields
        self.name = self.id #alias

    def __str__(self):
        return self.id

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

    def export(self):
        s = self.chrom + '\t' + str(self.begin) + '\t' + str(self.end)
        for k, v in self.fields.items():
            s += '\t' + getattr(self, k)
        s += '\n'
        return s

    def bedify(self):
        s = str(self.chrom) + '\t' + str(self.begin) + '\t' + str(self.end)
        for (k, v) in self.fields:
            s += '\t' + str(v)
        s += '\n'
        return s

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

def track_from_name(name):
    tokens = name.split('_')
    return BedTrack(tokens[0], int(tokens[1]), int(tokens[2]))

# keywords is an array of tuples: (name, default, transformation)
def load_tracks_from_file(path, parse_header = False, keywords = []):
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

def load_tracks_from_file_as_dict(path, parse_header = False, keywords = []):
    return {str(track): track for track in load_tracks_from_file(path, parse_header, keywords)}


