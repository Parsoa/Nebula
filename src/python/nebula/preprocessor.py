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
    vcf,
    config,
    counter,
    junction,
    reduction,
    map_reduce,
    statistics,
    visualizer,
    programming
)

from nebula.debug import *
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
        user_print('Loaded', len(tracks), 'tracks.')
        tracks = bed.sort_tracks(tracks)
        svtypes = {track.svtype: True for track in tracks}
        merged_tracks = []
        for svtype in svtypes:
            if svtype in TrackPreprocessorJob.SUPPORTED_SVTYPES:
                _tracks = [track for track in tracks if track.svtype == svtype]
                print('Merging', len(_tracks), svtype, 'tracks..')
                _tracks = bed.filter_overlapping_tracks(_tracks, svtype)
                print('Done', len(_tracks), 'tracks remaining.')
                merged_tracks += _tracks
        tracks = bed.sort_tracks(bed.filter_short_tracks(merged_tracks))
        tracks = {track.id: track for track in tracks if track.svtype in TrackPreprocessorJob.SUPPORTED_SVTYPES}
        user_print('Removed overlapping tracks.', len(tracks), 'non-overlapping tracks remainig.')
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

    def plot(self):
        values = open(os.path.join(self.get_current_job_directory(), 'depth_all.bed'), 'r').readlines()
        values = [float(line) for line in values]
        visualizer.scatter(list(range(0, 100)), values, 'GC_content_bp', self.get_current_job_directory(), 'GC', 'Coverage') 
        exit()

    def load_inputs(self):
        c = config.Configuration()
        self.gc = {}
        self.window_size = 500
        self.bin_size = self.window_size // 100
        for i in range(0, 100 + 1):
            self.gc[i] = {}
        self.load_reference_counts_provider() 
        self.chromosomes = extract_whole_genome()
        self.load_exon_regions()
        self.round_robin(self.chromosomes)

    def load_exon_regions(self):
        c = config.Configuration()
        exons = bed.load_tracks_from_file_as_dict(c.bed[0])
        self.exons = {}
        for chrom in self.chromosomes:
            self.exons[chrom] = {}
        for exon in exons:
            chrom = exons[exon].chrom.lower()
            if chrom in self.exons:
                self.exons[chrom][exon] = exons[exon]
        print(len(self.exons))

    # No exon, unique kmers from all across the genome
    def transform(self, sequence, chrom):
        c = config.Configuration()
        telomere_length = 10000
        sequence = sequence[telomere_length: -telomere_length]
        t = time.time()
        l = len(sequence)
        i = self.window_size // 2
        while i < l - self.window_size // 2:
            kmer = canonicalize(sequence[i - c.ksize // 2: i + c.ksize // 2])
            if self.reference_counts_provider.get_kmer_count(kmer) == 1:
                gc = calculate_gc_content(sequence[i - self.window_size // 2: i + self.window_size // 2]) // self.bin_size
                contig = chrom + ':' + str(telomere_length + i - self.window_size // 2) + '-' + str(telomere_length + i + self.window_size // 2)
                self.gc[gc][kmer] = {'contig': contig, 'gc': gc, 'coverage': gc}
                i += self.window_size // 2 + c.ksize // 2
            else:
                i += 1
                continue
            if i % 10000 == 0:
                s = time.time()
                p = i / float(len(sequence))
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
        gc_files = []
        for batch in self.load_output():
            for gc in range(0, 100 + 1):
                for kmer in batch[str(gc)]:
                    self.gc[gc][kmer] = batch[str(gc)][kmer]
        for gc in self.gc:
            print(len(self.gc[gc]), 'kmers with GC content', gc)
            if len(self.gc[gc]) > 10000:
                k = list(self.gc[gc].keys())
                for kmer in [k[i] for i in sorted(random.sample(xrange(len(k)), 10000))]:
                    self.kmers[kmer] = self.gc[gc][kmer]
            else:
                for kmer in self.gc[gc]:
                    self.kmers[kmer] = self.gc[gc][kmer]
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
        system_print_high('Counting', green(len(self.inner_kmers)), 'inner kmers')
        system_print_high('Counting', green(len(self.junction_kmers)), 'junction kmers')
        system_print_high('Counting', green(len(self.depth_kmers)), 'depth-of-coverage kmers')
        system_print_high('Counting', green(len(self.gc_kmers)), 'GC-content kmers')

    def merge_kmers(self):
        kmers = {}
        kmers['gc_kmers'] = self.gc_kmers
        kmers['depth_kmers'] = self.depth_kmers
        kmers['inner_kmers'] = self.inner_kmers
        kmers['junction_kmers'] = self.junction_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

    def load_inner_kmers(self):
        c = config.Configuration()
        self.inner_kmers = {}
        job = reduction.FilterInnerKmersJob()
        with open(os.path.join(job.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
            self.inner_kmers.update(json.load(json_file))

    def load_junction_kmers(self):
        c = config.Configuration()
        job = junction.FilterJunctionKmersJob()
        self.junction_kmers = {}
        name = 'kmers_cpp.json' if c.cpp else 'kmers.json'
        with open(os.path.join(job.get_current_job_directory(), name), 'r') as json_file:
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
        with open(os.path.join(self.get_current_job_directory(), 'reference_kmers.json'), 'w') as json_file:
            json.dump(self.depth_kmers, json_file, indent = 4)

    def load_gc_content_kmers(self):
        c = config.Configuration()
        job = GcContentKmerSelectionJob()
        gc_kmers = json.load(open(c.gckmers)) if c.gckmers else job.execute()
        self.gc_kmers = {}
        for kmer in gc_kmers:
            if kmer not in self.inner_kmers and kmer not in self.junction_kmers:
                self.gc_kmers[kmer] = gc_kmers[kmer]
                self.gc_kmers[kmer].update({'loci': {}, 'count': 0})

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
        system_print_high('Kmers exported for', len(self.tracks), 'tracks.')
        with open(os.path.join(self.get_current_job_directory(), 'tracks.bed'), 'w') as bed_file:
            for track in self.tracks:
                bed_file.write(c.tracks[track].serialize())
