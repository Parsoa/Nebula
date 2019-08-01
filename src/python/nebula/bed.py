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
        for (k, v) in fields:
            setattr(self, k, v)
        self.id = self.svtype + '@' + self.chrom + '_' + str(self.begin) + '_' + str(self.end)
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

    @staticmethod
    def json_serialize(self):
        d = {}
        d['chrom'] = self.chrom
        d['begin'] = self.begin
        d['end'] = self.end
        for (k, v) in self.fields:
            d[k] = v
        return d

    def add_field(self, key, value):
        setattr(self, key, value)
        self.fields.append((key, value))

    def serialize(self):
        s = ''
        s += str(self.chrom) + '\t' + str(self.begin) + '\t' + str(self.end)
        for (k, v) in self.fields:
            s += '\t' + v
        s += '\n'
        return s

    def header(self):
        s = '#CHROM\tBEGIN\tEND'
        for (key, value) in self.fields:
            s += '\t' + key
        s += '\n'
        return s

    def __getitem__(self, key):
        return getattr(self, key)

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

def track_from_name(name):
    tokens = name.split('_')
    return BedTrack(tokens[0], int(tokens[1]), int(tokens[2]))

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


