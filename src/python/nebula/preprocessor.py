from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import argparse
import operator
import traceback
import subprocess

from nebula import (
    bed,
    config,
    gapped,
    counter,
    junction,
    reduction,
    map_reduce,
    statistics,
    visualizer,
    programming
)

from nebula.kmers import *
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class TrackPreprocessorJob(map_reduce.Job):

    _name = 'TrackPreprocessorJob'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = TrackPreprocessorJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def execute(self):
        c = config.Configuration()
        tracks = {}
        for path in c.bed:
            tracks.update(bed.load_tracks_from_file_as_dict(path, parse_header = True))
        tracks = [tracks[track] for track in tracks]
        print('loaded', len(tracks), 'tracks')
        #tracks = self.filter_overlapping_tracks(\
        #            sorted(sorted(tracks, key = lambda x: x.begin), key = lambda y: y.chrom)\
        #        )
        tracks = {track.id: track for track in tracks}
        print('removed overlapping tracks:', len(tracks))
        return tracks

    def filter_overlapping_tracks(self, tracks):
        remove = []
        i = 0
        while i < len(tracks):
            for j in range(i + 1, len(tracks)):
                # j is contained inside i
                if tracks[j].chrom != tracks[i].chrom:
                    i = j
                    break
                if tracks[j].begin <= tracks[i].end:
                    remove.append(j)
                    print(red(str(tracks[j])), 'overlaps', blue(str(tracks[i])))
                    continue
                else:
                    i = j
                    break
            if i == len(tracks) - 1:
                break
        n = 0
        for index in sorted(remove):
            tracks.pop(index - n)
            n = n + 1
        return tracks

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class MixKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'MixKmersJob'
    _category = 'preprocessing'
    _previous_job = None 
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = MixKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.load_kmers()
        self.merge_kmers()
        self.export_tracks()
        exit()

    def load_kmers(self):
        c = config.Configuration()
        self.load_inner_kmers()
        self.load_junction_kmers()
        self.load_gc_content_kmers()
        self.load_depth_of_coverage_kmers()
        print('Counting', green(len(self.inner_kmers)), 'inner kmers')
        print('Counting', green(len(self.junction_kmers)), 'junction kmers')
        print('Counting', green(len(self.depth_kmers)), 'depth-of-coverage kmers')
        print('Counting', green(len(self.gc_kmers)), 'GC-content kmers')

    def merge_kmers(self):
        kmers = {}
        kmers['gc_kmers'] = self.gc_kmers
        kmers['depth_kmers'] = self.depth_kmers
        kmers['inner_kmers'] = self.inner_kmers
        kmers['junction_kmers'] = self.junction_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

    def load_inner_kmers(self):
        self.inner_kmers = {}
        job = reduction.FilterLociIndicatorKmersJob()
        with open(os.path.join(job.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
            self.inner_kmers.update(json.load(json_file))

    def load_junction_kmers(self):
        c = config.Configuration()
        job = junction.FilterJunctionKmersJob()
        self.junction_kmers = {}
        with open(os.path.join(job.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                kmer = canonicalize(kmer)
                if kmer in self.inner_kmers:
                    self.inner_kmers.pop(kmer, None)
                else:
                    self.junction_kmers[kmer] = kmers[kmer]

    def load_depth_of_coverage_kmers(self):
        n = 100000
        self.depth_kmers = {}
        self.load_reference_counts_provider() 
        for kmer, count in self.reference_counts_provider.stream_kmers():
            if count == 1 and kmer.find('N') == -1:
                canon = canonicalize(kmer)
                if not canon in self.inner_kmers and not canon in self.junction_kmers:
                    self.depth_kmers[canon] = {'loci': {}, 'count': 0}
                    n -= 1
                    if n == 0:
                        break
        self.unload_reference_counts_provider()

    def load_gc_content_kmers(self):
        self.gc_kmers = {}
    #    job = depth.ChromosomeGcContentEstimationJob()
    #    kmers = job.execute()
    #    for kmer in kmers:
    #        if kmer not in self.junction_kmers and kmer not in self.inner_kmers:
    #            self.gc_kmers[kmer] = {'gc': {}, 'loci': {}, 'count': 0}

    def export_tracks(self):
        c = config.Configuration()
        self.tracks = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.junction_kmers:
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        print('Kmers exported for', len(self.tracks), 'tracks')
        #for track in self.tracks:
        #    with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
        #        json.dump(self.tracks[track], json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'tracks.json'), 'w') as json_file:
        #    json.dump(self.tracks, json_file, indent = 4)
