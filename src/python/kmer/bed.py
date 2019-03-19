from __future__ import print_function

import io
import os
import pwd
import sys
import json
import time
import subprocess

from kmer import (
    config,
)

from kmer.kmers import *
from kmer.debug import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

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

    # TODO: This is really ugly, hopefully won't happen in practice
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
        for k, v in self.kwargs.items():
            s += '\t' + getattr(self, k)
        s += '\n'
        return s

    def extract_base_sequence(self):
        c = config.Configuration()
        # this is the largest sequence that we will ever need for this track
        # <- Slack -><-actual sequence-><- Slack ->
        self.slack = c.readlength - c.ksize
        begin = self.begin - self.slack
        end = self.end + self.slack
        chromosome = extract_chromosome(self.chrom.lower())
        self.sequence = chromosome[begin: end]
        if begin > len(chromosome) or end > len(chromosome) or len(self.sequence) == 0:
            print(self.chrom)
            print(len(chromosome))
            print(self.slack)
            print(self.begin)
            print(self.end)
            debug_breakpoint()

    def extract_inner_kmers(self, counter, count, n, overlap = True, canonical = True):
        c = config.Configuration()
        self.extract_base_sequence()
        begin = self.slack
        end = len(self.sequence) - self.slack
        inner_seq = self.sequence[begin : end + 1 - c.ksize]
        #print(inner_seq)
        #print(len(inner_seq))
        inner_kmers = c_extract_kmers(c.ksize, counter, count, overlap, canonical, inner_seq)
        if len(inner_kmers) <= n:
            return inner_kmers
        else:
            items = sorted(inner_kmers.items(), key = lambda item: item[1])[0:n]
            return {item[0]: item[1] for item in items}

    def extract_boundary_gapped_kmers(self, counter = lambda x: 1, count = 1):
        c = config.Configuration()
        self.extract_base_sequence()
        begin = self.slack
        end = len(self.sequence) - self.slack
        gapped_kmers = {}
        h = c.hsize
        #
        b = self.sequence[begin - h - 2 - c.ksize: begin + 3 + h + c.ksize]
        kmer = self.sequence[begin - h - 2: begin + 3 + h]
        prefix = self.sequence[begin - h - 2 - c.ksize: begin - h - 2]
        suffix = self.sequence[begin + 3 + h: begin + 3 + h + c.ksize]
        gapped_kmers[kmer] = {'left': prefix, 'right': suffix, 'side': 'inner'}
        #
        e = self.sequence[end - h - 2 - c.ksize: end + 3 + h + c.ksize]
        kmer = self.sequence[end - h - 2: end + 3 + h]
        prefix = self.sequence[end - h - 2 - c.size: end - h - 2]
        suffix = self.sequence[end + 3 + h: end + 3 + h + c.ksize]
        gapped_kmers[kmer] = {'left': prefix, 'right': suffix, 'side': 'inner'}
        #
        kmer = self.sequence[begin - 2 - h: begin + 3] + self.sequence[end - 2: end + 3 + h]
        prefix = self.sequence[begin - h - 2 - c.ksize: begin - h - 2]
        suffix = self.sequence[end + 3 + h: end + 3 + h + c.ksize]
        gapped_kmers[kmer] = {'left': prefix, 'right': suffix, 'side': 'outer'}
        return gapped_kmers

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


