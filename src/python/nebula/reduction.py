from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import time
import argparse
import operator
import traceback
import subprocess

from nebula import (
    bed,
    depth,
    config,
    counter,
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
)
from nebula.debug import *
from nebula.kmers import *
from nebula.commons import *
from nebula.chromosomes import *

import pysam
import numpy as np

print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractInnerKmersJob(map_reduce.Job):

    _name = 'ExtractInnerKmersJob'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractInnerKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def execute(self):
        c = config.Configuration()
        self.pre_process()
        self.create_output_directories()
        self.find_thread_count()
        self.load_inputs()

    def load_inputs(self):
        c = config.Configuration()
        extract_whole_genome()
        self.tracks = c.tracks
        self.load_reference_counts_provider()
        self.contigs = pysam.AlignmentFile(c.contigs, "rb") if c.contigs else None
        if self.contigs:
            self.extract_kmers()
        else:
            print('No assembly contigs provided, extracting from reference')
            self.round_robin(self.tracks, filter_func = lambda track: track.end - track.begin > 1000000)
            self.distribute_workload()
            self.wait_for_children()
            self.reduce()

    def extract_kmers(self):
        c = config.Configuration()
        contig_index = {}
        output = {}
        for track in self.tracks:
            contig_index[self.tracks[track].contig] = self.tracks[track]
        n = 0
        print(len(contig_index), 'total tracks')
        for read in self.contigs.fetch():
            if read.query_name in contig_index:
                track = contig_index[read.query_name]
                inner_kmers = self.extract_inner_kmers(track, read.query_sequence[int(track.contig_start) - c.ksize: int(track.contig_end) + c.ksize])
                if inner_kmers:
                    path = os.path.join(self.get_current_job_directory(), track.id + '.json')
                    with open(path, 'w') as json_file:
                        json.dump(inner_kmers, json_file, indent = 4)
                    output[track.id] = path
                n += 1
                if n % 1000 == 0:
                    print(n)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump(output, json_file, indent = 4)

    def transform(self, track, track_name):
        debug_breakpoint()
        print(green(track_name))
        c = config.Configuration()
        inner_kmers = self.extract_inner_kmers(track, None)
        if len(inner_kmers) == 0:
            print(red('skipping', track_name, 'no inner kmers found'))
            return None
        name = track_name  + '.json'
        with open(os.path.join(self.get_current_job_directory(), name), 'w') as json_file:
            json.dump(inner_kmers, json_file, sort_keys = True, indent = 4)
        return name

    def extract_inner_kmers(self, track, assembly):
        c = config.Configuration()
        chrom = extract_chromosome(track.chrom.lower())
        if not chrom:
            print(red(track.chrom.lower()))
            return None
        inner_kmers = {}
        if track.svtype == 'INS':
            if assembly:
                inner_kmers = self.extract_insertion_kmers(track, assembly)
            elif hasattr(track, 'seq'):
                inner_kmers = self.extract_insertion_kmers(track, track.seq.upper())
        if track.svtype == 'DEL':
            inner_kmers = self.extract_deletion_kmers(track, chrom)
        # No inner kmers for INV, MEI or ALU
        if len(inner_kmers) > 1000:
            items = sorted(inner_kmers.items(), key = lambda item: item[1]['reference'])[0 : 1000]
            return {item[0]: item[1] for item in items}
        return inner_kmers

    def extract_insertion_kmers(self, track, inner_seq):
        c = config.Configuration()
        inner_kmers = {}
        if len(inner_seq) >= 3 * c.ksize:
            for i in range(c.ksize, len(inner_seq) - 2 * c.ksize + 1):
                kmer = canonicalize(inner_seq[i: i + c.ksize])
                if not kmer in inner_kmers:
                    inner_kmers[kmer] = {
                        'loci': {},
                        'tracks': {},
                        'reference': self.reference_counts_provider.get_kmer_count(kmer),
                    }
                inner_kmers[kmer]['loci']['inside_' + track.id + '@' + str(i)] = {
                    'seq': {
                        'all': inner_seq[i - c.ksize: i + c.ksize + c.ksize],
                        'left': inner_seq[i - c.ksize: i],
                        'right': inner_seq[i + c.ksize: i + c.ksize + c.ksize]
                    }
                }
                if not track.id in inner_kmers[kmer]['tracks']:
                    inner_kmers[kmer]['tracks'][track.id] = 0
                inner_kmers[kmer]['tracks'][track.id] += 1
        return inner_kmers

    def extract_deletion_kmers(self, track, chrom):
        c = config.Configuration()
        inner_kmers = {}
        slack = 5 # arbitray
        inner_seq = chrom[track.begin + slack: track.end - slack]
        for i in range(len(inner_seq) - c.ksize + 1):
            kmer = canonicalize(inner_seq[i: i + c.ksize])
            if not kmer in inner_kmers:
                inner_kmers[kmer] = {
                    'loci': {},
                    'tracks': {},
                    'reference': self.reference_counts_provider.get_kmer_count(kmer),
                }
            if not track.id in inner_kmers[kmer]['tracks']:
                inner_kmers[kmer]['tracks'][track.id] = 0
            inner_kmers[kmer]['tracks'][track.id] += 1
        return inner_kmers

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractLociIndicatorKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'ExtractLociIndicatorKmersJob'
    _category = 'preprocessing'
    _previous_job = ExtractInnerKmersJob

    @staticmethod
    def launch(**kwargs):
        job = ExtractLociIndicatorKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.chroms = extract_whole_genome()
        self.tracks = self.load_previous_job_results()
        self.inner_kmers = {}
        for track in self.tracks:
            with open(os.path.join(self.get_previous_job_directory(), self.tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                self.keep_best_kmers(kmers)
        print('Finding loci for', green(len(self.inner_kmers)), 'kmers')
        self.round_robin(self.chroms)

    # This will lose some shared kmers
    def keep_best_kmers(self, kmers):
        k = {}
        for i in range(0, 5 + 1):
            k[i] = []
        for kmer in kmers:
            ref = kmers[kmer]['reference']
            if ref <= 5:
                k[ref].append(kmer)
        n = 0
        for i in range(0, 5 + 1):
            for kmer in k[i]:
                if n == 50:
                    break
                if not kmer in self.inner_kmers:
                    self.inner_kmers[kmer] = kmers[kmer]
                self.inner_kmers[kmer]['loci'].update(kmers[kmer]['loci'])
                self.inner_kmers[kmer]['tracks'].update(kmers[kmer]['tracks'])
                n += 1

    def transform(self, sequence, chrom):
        c = config.Configuration()
        index = 0
        slack = c.ksize
        t = time.time()
        for kmer in stream_kmers(c.ksize, True, True, sequence):
            if kmer in self.inner_kmers:
                locus = chrom + '_' + str(index)
                self.inner_kmers[kmer]['loci'][locus] = {
                    'seq': {
                        'all': sequence[index - slack: index + c.ksize + slack],
                        'left': sequence[index - slack: index],
                        'right': sequence[index + c.ksize: index + c.ksize + slack]
                    }
                }
            index += 1
            if index % 10000 == 0:
                s = time.time()
                p = index / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
        return None

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.inner_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            if not os.path.isfile(path):
                print(red('couldn\'t find batch'), i, red('results will be unreliable'))
                continue
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    self.inner_kmers[kmer]['loci'].update(batch[kmer]['loci'])
        for track in self.tracks:
            self.tracks[track] = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                self.tracks[track][kmer] = self.inner_kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'w') as track_file:
                json.dump(self.tracks[track], track_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'indicator_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

    def plot_kmer_reference_count(self):
        x = []
        for track in self.tracks:
            print(track)
            with open(os.path.join(self.get_current_job_directory(), self.tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)['inner_kmers']
                for kmer in kmers:
                    x.append(kmers[kmer]['reference'])
        visualizer.histogram(x = x, name = 'reference_count', path = self.get_current_job_directory(), x_label = 'kmer count in reference', y_label = 'number of kmers', step = 0.1)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class FilterLociIndicatorKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'FilterLociIndicatorKmersJob'
    _category = 'preprocessing'
    _previous_job = ExtractLociIndicatorKmersJob

    @staticmethod
    def launch(**kwargs):
        job = FilterLociIndicatorKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = self.load_previous_job_results()
        self.round_robin(self.tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            t = c.tracks[track_name]
            for kmer in kmers:
                if kmer.find('N') != -1:
                    continue
                seed = sum(list(map(lambda s: ord(s), kmer)))
                if not kmer in self.kmers:
                    self.kmers[kmer] = {}
                    self.kmers[kmer]['loci'] = kmers[kmer]['loci']
                    self.kmers[kmer]['total'] = 0
                    self.kmers[kmer]['count'] = 0
                    self.kmers[kmer]['doubt'] = 0
                    self.kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                    self.kmers[kmer]['reference'] = kmers[kmer]['reference']
                    self.kmers[kmer]['interest_masks'] = {}
                    for locus in self.kmers[kmer]['loci']:
                        self.kmers[kmer]['loci'][locus]['masks'] = {self.kmers[kmer]['loci'][locus]['seq']['left']: True, self.kmers[kmer]['loci'][locus]['seq']['right']: True}
                for locus in self.kmers[kmer]['loci']:
                    tokens = locus.split('_')
                    if 'inside' in locus or (tokens[0].lower() == t.chrom.lower() and int(tokens[1]) >= t.begin and int(tokens[1]) < t.end):
                        self.kmers[kmer]['interest_masks'].update(self.kmers[kmer]['loci'][locus]['masks'])
        return None

    def is_kmer_returning(self, kmer):
        c = config.Configuration()
        for track in kmer['tracks']:
            t = c.tracks[track]
            for loci in kmer['loci']:
                if not 'inside' in loci:
                    start = int(loci.split('_')[1])
                    if abs(start - t.begin) < 2 * c.ksize:
                        return True
                    if abs(start - t.end) < 2 * c.ksize:
                        return True
        return False

    def output_batch(self, batch):
        for kmer in self.kmers:
            for locus in self.kmers[kmer]['loci'].keys():
                l = self.get_shared_masks(self.kmers[kmer]['interest_masks'], self.kmers[kmer]['loci'][locus]['masks'])
                if l == 0:
                    self.kmers[kmer]['loci'].pop(locus, None)
            self.kmers[kmer].pop('interest_masks', None)
            self.kmers[kmer]['reduction'] = self.kmers[kmer]['reference']
            self.kmers[kmer]['reference'] = len(self.kmers[kmer]['loci'])
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def get_shared_masks(self, interest_kmers, kmers):
        l = []
        for x in kmers:
            for y in interest_kmers:
                if is_canonical_subsequence(x[4:-4], y):
                    l.append(x)
                    break
        return len(l)

    def reduce(self):
        self.kmers = {}
        for batch in self.load_output():
            for kmer in batch:
                self.kmers[kmer] = batch[kmer]
        #
        print(len(self.kmers), 'kmers')
        returning_kmers = {kmer: self.kmers[kmer] for kmer in self.kmers if self.is_kmer_returning(self.kmers[kmer])}
        with open(os.path.join(self.get_current_job_directory(), 'returning.json'), 'w') as json_file:
            json.dump(returning_kmers, json_file, indent = 4)
        self.kmers = {kmer: self.kmers[kmer] for kmer in self.kmers if not self.is_kmer_returning(self.kmers[kmer])}
        print(len(self.kmers), 'kmers after filtering returning kmers')
        #
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        #
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class AdjustCoverageGcContentJob(map_reduce.BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'AdjustCoverageGcContentJob'
    _category = 'preprocessing'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = AdjustCoverageGcContentJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            self.kmers = json.load(json_file)
        print(len(self.kmers))
        gc_coverage_job = depth.ChromosomeGcContentEstimationJob()
        with open(os.path.join(gc_coverage_job.get_current_job_directory(), 'coverage.json'), 'r') as json_file:
            self.coverage = json.load(json_file)
        print('Adjusting GC coverage for', green(len(self.kmers)), 'kmers')
        n = 0
        for kmer in self.kmers:
            self.transform(self.kmers[kmer], kmer)
            n += 1
            if n % 1000 == 0:
                print(n, 'out of', len(self.kmers))
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        for track in self.tracks:
            print('exporting track', track)
            with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'indicator_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)
        exit()

    def transform(self, kmer, k):
        coverage = []
        for locus in kmer['loci']:
            gc = calculate_gc_content(kmer['loci'][locus]['seq']['all'])
            kmer['loci'][locus]['coverage'] = self.coverage[gc]
            coverage.append(self.coverage[gc])
        kmer['coverage'] = statistics.mean(coverage)

