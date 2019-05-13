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
        self.half_mers = {}
        self.depth_kmers = {}
        self.inner_kmers = {}
        self.gapped_kmers = {}
        self.junction_kmers = {}
        self.load_kmers()
        self.merge_kmers()
        #self.export_tracks()
        exit()

    def load_kmers(self):
        c = config.Configuration()
        #job = reduction.FilterLociIndicatorKmersJob()
        self.inner_kmers = {}
        #with open(os.path.join(job.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
        #    self.inner_kmers.update(json.load(json_file))
        self.load_junction_kmers()
        self.load_depth_of_coverage_kmers()
        #self.load_gapped_kmers()

    def merge_kmers(self):
        kmers = {}
        kmers['half_mers'] = self.half_mers
        kmers['depth_kmers'] = self.depth_kmers
        kmers['inner_kmers'] = self.inner_kmers
        kmers['gapped_kmers'] = self.gapped_kmers
        kmers['junction_kmers'] = self.junction_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

    def load_junction_kmers(self):
        c = config.Configuration()
        job = junction.CountJunctionKmersJob()
        with open(os.path.join(job.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                kmer = canonicalize(kmer)
                if kmer in self.inner_kmers:
                    self.inner_kmers.pop(kmer, None)
                else:
                    self.junction_kmers[kmer] = kmers[kmer]
        #with open(os.path.join(self.get_current_job_directory(), 'junction_kmers.json'), 'w') as json_file:
        #    json.dump(self.junction_kmers, json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
        #    json.dump(self.inner_kmers, json_file, indent = 4)

    def load_depth_of_coverage_kmers(self):
        n = 100000
        self.load_reference_counts_provider() 
        for kmer, count in self.reference_counts_provider.stream_kmers():
            if count == 1 and kmer.find('N') == -1:
                canon = canonicalize(kmer)
                if not canon in self.inner_kmers and not canon in self.junction_kmers:
                    self.depth_kmers[canon] = {'loci': {}, 'count': 0}
                    n -= 1
                    if n == 0:
                        break
        print('Counting', green(len(self.depth_kmers)), 'depth kmers')
        #with open(os.path.join(self.get_current_job_directory(), 'depth_kmers.json'), 'w') as json_file:
        #    json.dump(self.depth_kmers, json_file, sort_keys = True, indent = 4)
        self.unload_reference_counts_provider()

    def load_gapped_kmers(self):
        c = config.Configuration()
        job = gapped.CountUniqueGappedKmersJob()
        tracks = job.load_previous_job_results()
        for track in tracks:
            with open(os.path.join(job.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers:
                    if kmers[kmer]['gap'] != -1:
                        left = kmer[:c.hsize]
                        right = kmer[-c.hsize:]
                        self.gapped_kmers[kmer] = kmers[kmer]
                        self.gapped_kmers[kmer]['count'] = 0
                        self.gapped_kmers[kmer]['doubt'] = 0
                        if not left in self.half_mers:
                            self.half_mers[left] = {}
                        self.half_mers[left][right] = kmer 
                        left = reverse_complement(left)
                        right = reverse_complement(right)
                        if not right in self.half_mers:
                            self.half_mers[right] = {}
                        self.half_mers[right][left] = kmer
        print('Counting', green(len(self.gapped_kmers)), 'gapped kmers')
        #with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
        #    json.dump(self.half_mers, json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
        #    json.dump(self.gapped_kmers, json_file, indent = 4)

    def export_tracks(self):
        c = config.Configuration()
        self.tracks = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.gapped_kmers:
            for track in self.gapped_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['gapped_kmers'][kmer] = self.gapped_kmers[kmer]
                self.tracks[track]['gapped_kmers'][kmer]['type'] = 'gapped'
        for kmer in self.junction_kmers:
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'tracks.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, indent = 4)
