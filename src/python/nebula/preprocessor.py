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
    counter,
    junction,
    reduction,
    map_reduce,
    statistics,
    visualizer,
    programming
)

from nebula.kmers import *
from nebula.logger import *
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

    SUPPORTED_SVTYPES = ['DEL', 'INS', 'INV', 'ALU', 'MEI']

    def execute(self):
        c = config.Configuration()
        tracks = {}
        for path in c.bed:
            tracks.update(bed.load_tracks_from_file_as_dict(path, parse_header = True))
        tracks = [tracks[track] for track in tracks]
        user_print('Loaded', len(tracks), 'tracks.')
        if self.filter_overlap:
            tracks = self.filter_overlapping_tracks(\
                        sorted(sorted(tracks, key = lambda x: x.begin), key = lambda y: y.chrom)\
                    )
        tracks = {track.id: track for track in tracks if track.svtype in TrackPreprocessorJob.SUPPORTED_SVTYPES}
        user_print('Removed overlapping tracks.', len(tracks), ' non-overlapping tracks.')
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
                    user_print_warning(str(tracks[j]), 'overlaps', blue(str(tracks[i])))
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

class GcContentKmerSelectionJob(map_reduce.GenomeDependentJob):

    _name = 'ChromosomeGcContentEstimationJob'
    _category = 'preprocessing'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = ChromosomeGcContentEstimationJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.gc = {}
        self.window_size = 96
        for i in range(0, self.window_size + 1):
            self.gc[i] = {}
        self.chromosomes = extract_whole_genome()
        self.load_reference_counts_provider() 
        self.round_robin(self.chromosomes)

    def transform(self, sequence, chrom):
        l = len(sequence)
        d = l / self.window_size
        t = time.time()
        for i in range(0, d):
            window = sequence[i * self.window_size: (i + 1) * self.window_size]
            gc = calculate_gc_content(window)
            kmers = c_extract_kmers(32, self.reference_counts_provider.get_kmer_count, 1, True, True, window)
            for kmer in kmers:
                self.gc[gc][kmer] = 0
                break
            if i % 100 == 0:
                s = time.time()
                p = (len(sequence) - i * 100) / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
        return None

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.gc, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def reduce(self):
        c = config.Configuration()
        self.kmers = {}
        random.seed(c.seed)
        for batch in self.load_output():
            for i in range(0, self.window_size + 1):
                for kmer in batch[str(i)]:
                    self.gc[i][kmer] = 0
        for gc in self.gc:
            print(len(self.gc[gc]), 'kmers with GC content', gc)
            if len(self.gc[gc]) > 10000:
                k = list(self.gc[gc].keys())
                for kmer in [k[i] for i in sorted(random.sample(xrange(len(k)), 10000))]:
                    self.kmers[kmer] = gc
            else:
                for kmer in self.gc[gc]:
                    self.kmers[kmer] = gc
        print('counting', len(self.kmers), 'kmers')
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        return self.gc

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
        c = config.Configuration()
        job = GcContentKmerSelectionJob()
        gc_kmers = json.load(open(c.gckmers)) if c.gckmers else job.execute()
        self.gc_kmers = {}
        for kmer in gc_kmers:
            if kmer not in self.inner_kmers and kmer not in self.junction_kmers:
                self.gc_kmers[kmer] = {'loci': {}, 'gc': gc_kmers[kmer], 'count': 0}

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
        with open(os.path.join(self.get_current_job_directory(), 'tracks.bed'), 'w') as bed_file:
            for track in self.tracks:
                bed_file.write(c.tracks[track].serialize())
        #for track in self.tracks:
        #    with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
        #        json.dump(self.tracks[track], json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'tracks.json'), 'w') as json_file:
        #    json.dump(self.tracks, json_file, indent = 4)
