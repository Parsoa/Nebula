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

class VcfTrack:

    def __init__(self, preamble, fields = []):
        self.preamble = preamble
        self.fields = []
        self.genotype = '0/0'
        for (k, v) in fields:
            self[k] = v
            self.fields.append(k)
        info = self.info
        tokens = info.split(';')
        for token in tokens:
            if '=' in token:
                name, value = token.split('=')
                self[name] = value
        self.begin = int(self.pos)

    def serialize(self):
        s = ''
        for key in self.fields:
            s += str(self[key]) + '\t'
        s = s.strip()
        s += '\n'
        return s

    def header(self):
        s = ''
        for p in self.preamble:
            s += p + '\n'
        s += '#'
        for key in self.fields:
            s += key.upper() + '\t'
        s = s.strip()
        s += '\n'
        return s

    def __str__(self):
        return self.chrom + '_' + self.pos + '_' + self.end

    def __getitem__(self, key):
        key = key.lower()
        return getattr(self, key)

    def __setitem__(self, key, value):
        key = key.lower()
        setattr(self, key, value)

    def __delattr__(self, key):
        key = key.lower()
        del self.__dict__[key]

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

# keywords is an array of tuples: (name, default, transformation)
def load_tracks_from_file(path):
    tracks = []
    header = ''
    preamble = []
    with open(path, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('##'):
                preamble.append(line.strip())
            elif line.startswith('#'):
                header = line
            else:
                tokens = line.split()
                names = header[1:].split()
                fields = []
                for i, t in enumerate(tokens):
                    fields.append((names[i] if i != 9 else 'genotype', t))
                tracks.append(VcfTrack(preamble, fields))
            line = f.readline()
    return tracks

def load_tracks_from_file_as_dict(path, parse_header = True, keywords = []):
    return {str(track): track for track in load_tracks_from_file(path)}


