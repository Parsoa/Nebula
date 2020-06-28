from __future__ import print_function

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
from nebula.logger import *
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

    def load_inputs(self):
        c = config.Configuration()
        self.chroms = extract_whole_genome()
        self.tracks = c.tracks
        self.load_reference_counts_provider()
        self.round_robin(self.tracks, filter_func = lambda track: track.end - track.begin > 1000000)

    def transform(self, track, track_name):
        debug_breakpoint()
        c = config.Configuration()
        # No inner kmers for INV, MEI or ALU
        if track.svtype == 'INV':
            system_print_warning('Skipping inner kmer extraction for inversion ' + track_name + '.')
            return None
        inner_kmers = self.extract_inner_kmers(track)
        if len(inner_kmers) == 0:
            system_print_warning('No inner kmers found for ' + track_name + '.')
            return None
        name = track_name  + '.json'
        with open(os.path.join(self.get_current_job_directory(), name), 'w') as json_file:
            json.dump(inner_kmers, json_file, sort_keys = True, indent = 4)
        return name

    def extract_inner_kmers(self, track):
        c = config.Configuration()
        chrom = self.chroms[track.chrom]
        if not chrom:
            system_print_error('Chromosome', track.chrom.lower(), 'not found. Skipping track.')
            return {}
        inner_kmers = {}
        if track.svtype == 'INS':
            if hasattr(track, 'seq'):
                inner_kmers = self.extract_insertion_kmers(track, track.seq.upper(), chrom)
        if track.svtype == 'DEL':
            inner_kmers = self.extract_deletion_kmers(track, chrom)
        if len(inner_kmers) > 1000:
            items = sorted(inner_kmers.items(), key = lambda item: item[1]['reference'])[0 : 1000]
            return {item[0]: item[1] for item in items}
        return inner_kmers

    def extract_insertion_kmers(self, track, inner_seq, chrom):
        c = config.Configuration()
        inner_kmers = {}
        seq = chrom[track.begin - 250: track.begin] + inner_seq + chrom[track.end: track.end + 250]
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
                    },
                    'gc': calculate_gc_content(seq[250 + i + c.ksize // 2 - 250: 250 + i + c.ksize // 2 + 250]) // 5
                }
                if not track.id in inner_kmers[kmer]['tracks']:
                    inner_kmers[kmer]['tracks'][track.id] = 0
                inner_kmers[kmer]['tracks'][track.id] += 1
        return inner_kmers

    def extract_deletion_kmers(self, track, chrom):
        c = config.Configuration()
        inner_kmers = {}
        returning_kmers = {}
        new_seq = chrom[track.begin - c.ksize: track.begin] + chrom[track.end: track.end + c.ksize]
        for i in range(len(new_seq) - c.ksize + 1):
            kmer = canonicalize(new_seq[i: i + c.ksize])
            returning_kmers[kmer] = True
        inner_seq = chrom[track.begin - c.ksize: track.end + c.ksize]
        for i in range(len(inner_seq) - c.ksize + 1):
            kmer = canonicalize(inner_seq[i: i + c.ksize])
            if not kmer in returning_kmers:
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

class ScoreInnerKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'ScoreInnerKmersJob'
    _category = 'preprocessing'
    _previous_job = ExtractInnerKmersJob

    @staticmethod
    def launch(**kwargs):
        job = ScoreInnerKmersJob(**kwargs)
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
        if len(self.inner_kmers) == 0:
            return None
        l = float(len(sequence))
        index = 0
        t = time.time()
        for kmer in stream_kmers(c.ksize, True, True, sequence):
            if kmer in self.inner_kmers:
                locus = chrom + '_' + str(index)
                self.inner_kmers[kmer]['loci'][locus] = {
                    'seq': {
                        'all': sequence[index - c.ksize: index + c.ksize + c.ksize],
                        'left': sequence[index - c.ksize: index],
                        'right': sequence[index + c.ksize: index + c.ksize + c.ksize]
                    },
                    'gc': calculate_gc_content(sequence[index + c.ksize // 2 - 250: index + c.ksize // 2 + 250]) // 5
                }
            index += 1
            if index % 1000000 == 0:
                s = time.time()
                p = index / l
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
        return None

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.inner_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        no_mask_kmers = {}
        for batch in self.load_output():
            for kmer in batch:
                self.inner_kmers[kmer]['loci'].update(batch[kmer]['loci'])
        for track in self.tracks:
            self.tracks[track] = {}
        for kmer in self.inner_kmers:
            if len(self.inner_kmers[kmer]['loci']) == 0:
                no_mask_kmers[kmer] = self.inner_kmers[kmer]
            for track in self.inner_kmers[kmer]['tracks']:
                self.tracks[track][kmer] = self.inner_kmers[kmer]
        json.dump(no_mask_kmers, open(os.path.join(self.get_current_job_directory(), 'no_mask_kmers.json'), 'w'), indent = 4)
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as track_file:
                json.dump(self.tracks[track], track_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)
        return self.tracks

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class FilterInnerKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'FilterInnerKmersJob'
    _category = 'preprocessing'
    _previous_job = ScoreInnerKmersJob

    @staticmethod
    def launch(**kwargs):
        job = FilterInnerKmersJob(**kwargs)
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
            for kmer in kmers:
                if kmer.find('N') != -1:
                    continue
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
                    if 'inside_' in locus:
                        self.kmers[kmer]['interest_masks'].update(self.kmers[kmer]['loci'][locus]['masks'])
                        continue
                    for track in self.kmers[kmer]['tracks']:
                        t = c.tracks[track]
                        if tokens[0].lower() == t.chrom.lower() and int(tokens[1]) >= t.begin - c.ksize and int(tokens[1]) < t.end:
                            self.kmers[kmer]['interest_masks'].update(self.kmers[kmer]['loci'][locus]['masks'])
                            continue
        return None

    # Why only doing this for insertions?
    def is_kmer_returning(self, kmer):
        c = config.Configuration()
        for track in kmer['tracks']:
            t = c.tracks[track]
            if t.svtype == 'INS':
                for locus in kmer['loci']:
                    if not 'inside' in locus:
                        start = int(locus.split('_')[1])
                        if abs(start - t.end) < 2 * c.ksize:
                            return True
                        if abs(start - t.begin) < 2 * c.ksize:
                            return True
        return False

    def get_shared_masks(self, interest_kmers, kmers):
        l = []
        for x in kmers:
            for y in interest_kmers:
                if is_canonical_subsequence(x[4:-4], y):
                    l.append(x)
                    break
        return len(l)

    def output_batch(self, batch):
        # different batches may output the same kmer. Results will be the same.
        for kmer in self.kmers:
            self.kmers[kmer]['filtered_loci'] = {}
            for locus in list(self.kmers[kmer]['loci'].keys()):
                l = self.get_shared_masks(self.kmers[kmer]['interest_masks'], self.kmers[kmer]['loci'][locus]['masks'])
                if l == 0:
                    self.kmers[kmer]['filtered_loci'][locus] = self.kmers[kmer]['loci'][locus]
                    self.kmers[kmer]['loci'].pop(locus, None)
            self.kmers[kmer].pop('interest_masks', None)
            self.kmers[kmer]['reduction'] = self.kmers[kmer]['reference']
            self.kmers[kmer]['reference'] = len(self.kmers[kmer]['loci'])
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def reduce(self):
        c = config.Configuration()
        self.kmers = {}
        for batch in self.load_output():
            for kmer in batch:
                self.kmers[kmer] = batch[kmer]
        # make sure everything is valid
        for kmer in self.kmers:
            ts = [c.tracks[track] for track in self.kmers[kmer]['tracks']]
            for loci in self.kmers[kmer]['filtered_loci']:
                tokens = loci.split('_')
                begin = int(tokens[1])
                for t in ts:
                    assert tokens[0] != t.chrom or begin < t.begin - c.ksize or begin >= t.end
        # filter returning kmers
        system_print_high('Extracted', len(self.kmers), 'kmers.')
        returning_kmers = {kmer: self.kmers[kmer] for kmer in self.kmers if self.is_kmer_returning(self.kmers[kmer])}
        with open(os.path.join(self.get_current_job_directory(), 'returning.json'), 'w') as json_file:
            json.dump(returning_kmers, json_file, indent = 4)
        self.kmers = {kmer: self.kmers[kmer] for kmer in self.kmers if not self.is_kmer_returning(self.kmers[kmer])}
        system_print_high(len(self.kmers), 'kmers remaining after filtering returning kmers.')
        # dump all inner kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        # dump kmers for each track
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        system_print_high('Exported', len(self.tracks), 'tracks with eligible inner kmers.')
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)

