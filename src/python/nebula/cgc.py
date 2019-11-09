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

#from sets import Set

from nebula import (
    bed,
    config,
    counter,
    junction,
    reduction,
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
    preprocessor,
)

from nebula.debug import *
from nebula.kmers import *
from nebula.logger import *
from nebula.chromosomes import *

import numpy as np

from pulp import *
from scipy.signal import savgol_filter

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcCounterJob(counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'CgcCounterJob'
    _category = 'genotyping'
    _previous_job = preprocessor.MixKmersJob
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = CgcCounterJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.load_kmers()
        self.export_counter_input()
        self.round_robin()

    def load_kmers(self):
        c = config.Configuration()
        path = c.kmers[0]
        with open(path, 'r') as json_file:
            self.kmers = json.load(json_file)
        self.gc_kmers = self.kmers['gc_kmers']
        self.depth_kmers = self.kmers['depth_kmers']
        self.inner_kmers = self.kmers['inner_kmers']
        self.junction_kmers = self.kmers['junction_kmers']
        for path in c.kmers[1:]:
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                for kmer_type in ['junction_kmers', 'inner_kmers']:
                    x = getattr(self, kmer_type)
                    for kmer in kmers[kmer_type]:
                        if not kmer in x:
                            x[kmer] = kmers[kmer_type][kmer]
                        else:
                            for locus in kmers[kmer_type][kmer]['loci']:
                                if not locus in x[kmer]['loci']:
                                    x[kmer]['loci'][locus] = kmers[kmer_type][kmer]['loci']
                                else:
                                    x[kmer]['loci'][locus]['masks'].update(kmers[kmer_type][kmer]['loci'][locus]['masks'])
                            x[kmer]['tracks'].update(kmers[kmer_type][kmer]['tracks'])
        user_print('Counting:')
        user_print('Junction kmers:', blue(len(self.junction_kmers)))
        user_print('Inner kmers:', blue(len(self.inner_kmers)))
        user_print('Depth kmers:', blue(len(self.depth_kmers)))
        user_print('GC kmers:', blue(len(self.gc_kmers)))

    def export_counter_input(self):
        for kmer in self.gc_kmers:
            self.gc_kmers[kmer]['count'] = 0
        if self.resume_from_reduce:
            return
        _kmers = {}
        for kmers in [self.inner_kmers, self.junction_kmers]:
            for kmer in kmers:
                _kmers[kmer] = {}
                _kmers[kmer]['loci'] = {}
                _kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                for locus in kmers[kmer]['loci']:
                    _kmers[kmer]['loci'][locus] = {}
                    _kmers[kmer]['loci'][locus]['masks'] = kmers[kmer]['loci'][locus]['masks']
        for kmer in self.depth_kmers:
            _kmers[kmer] = {}
            _kmers[kmer]['loci'] = {}
        for kmer in self.gc_kmers:
            _kmers[kmer] = {}
            _kmers[kmer]['loci'] = {}
        system_print('Dumping', green(len(_kmers)), 'kmers for the counter')
        with open(os.path.join(self.get_current_job_directory(), 'pre_inner_kmers.json'), 'w') as json_file:
            json.dump(_kmers, json_file, indent = 4)

    def merge_count(self, kmer, tokens):
        count = tokens[0] 
        total = tokens[1] 
        canon = canonicalize(kmer)
        if canon in self.gc_kmers:
            self.gc_kmers[canon]['count'] += count / 2
        if canon in self.depth_kmers:
            self.depth_kmers[canon]['count'] += count / 2
        if canon in self.inner_kmers:
            self.inner_kmers[canon]['count'] += count / 2
            self.inner_kmers[canon]['total'] += total / 2
        if canon in self.junction_kmers:
            self.junction_kmers[canon]['count'] += count / 2
            self.junction_kmers[canon]['total'] += total / 2

    def estimate_depth_of_coverage(self):
        c = config.Configuration()
        counts = list(map(lambda kmer: self.depth_kmers[kmer]['count'], self.depth_kmers))
        visualizer.histogram(counts, 'depth_unfiltered', self.get_current_job_directory(), 'count', 'kmers', 1)
        print(len(counts), 'depth kmers')
        while True:
            std = np.std(counts)
            mean = np.mean(counts)
            print(mean, std)
            counts = list(filter(lambda x: x < 3 * mean, counts))
            if abs(np.mean(counts) - mean) < 0.5:
                break
        #
        visualizer.histogram(counts, 'depth_filtered', self.get_current_job_directory(), 'count', 'kmers', 1)
        user_print('Sample coverage estimate:', str(mean) + 'x')
        stats = {'coverage': int(round(mean)), 'std': std}
        with open(os.path.join(self.get_current_job_directory(), 'stats_' + str(c.ksize) + '.json'), 'w') as json_file:
            json.dump(stats, json_file, sort_keys = True, indent = 4)
        return stats

    def estimate_gc_content_coverage(self):
        c = config.Configuration()
        self.gc = [0] * (100 + 1)
        for gc in range(0, 100 + 1):
            counts = [self.gc_kmers[kmer]['count'] for kmer in self.gc_kmers if self.gc_kmers[kmer]['gc'] == gc]
            if len(counts) == 0 or not any(counts):
                self.gc[gc] = self.stats['coverage']
                continue
            while True:
                std = np.std(counts)
                mean = np.mean(counts)
                counts = list(filter(lambda x: x <= 3 * mean, counts))
                if abs(np.mean(counts) - mean) < 0.5:
                    self.gc[gc] = max(mean, 0.2 * self.stats['coverage'])
                    break
        visualizer.scatter(list(range(0, 100 + 1)), self.gc, 'GC_coverage', self.get_current_job_directory(), 'GC', 'coverage')
        self.gc = savgol_filter(self.gc, 7, 3)
        visualizer.scatter(list(range(0, 100 + 1)), self.gc, 'smooth_GC_coverage', self.get_current_job_directory(), 'GC', 'coverage')

    def adjust_kmer_gc_coverage(self, kmer, kmer_seq):
        c = config.Configuration()
        kmer['gc_coverage'] = sum(map(lambda locus: self.gc[kmer['loci'][locus]['gc']], kmer['loci'])) / len(kmer['loci'])

    def reduce(self):
        self.merge_counts()
        self.export_counted_kmers()
        self.stats = self.estimate_depth_of_coverage()
        self.estimate_gc_content_coverage()
        self.export_genotyping_tracks()
        return self.tracks, self.stats 

    def export_counted_kmers(self):
        c = config.Configuration()
        kmers = {}
        kmers['gc_kmers'] = self.gc_kmers
        kmers['depth_kmers'] = self.depth_kmers
        kmers['inner_kmers'] = self.inner_kmers
        kmers['junction_kmers'] = self.junction_kmers
        if c.cgc:
            name = 'kmers_' + (c.fastq.split('/')[-1] if c.fastq else c.bam.split('/')[-1]) + '.json'
            path = os.path.join(c.workdir, name)
        else:
            name = 'kmers.json'
            path = os.path.join(self.get_current_job_directory(), name)
        with open(path, 'w') as json_file:
            json.dump(kmers, json_file, indent = 4)

    def export_genotyping_tracks(self):
        c = config.Configuration()
        self.tracks = {}
        for kmer in self.inner_kmers:
            self.adjust_kmer_gc_coverage(self.inner_kmers[kmer], kmer)
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.junction_kmers:
            self.junction_kmers[kmer]['gc_coverage'] = self.stats['coverage']
            if 'inverse' in self.junction_kmers[kmer]:
                self.junction_kmers[kmer]['count'] = self.junction_kmers[kmer]['total'] - self.junction_kmers[kmer]['count']
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'CgcIntegerProgrammingJob'
    _category = 'genotyping'
    _previous_job = CgcCounterJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.round_robin(self.tracks)
        self.resume_from_reduce = False
        self.lp_kmers = {}

    def transform(self, kmers, track_name):
        c = config.Configuration()
        lp_kmers = {}
        all_non_unique = all(map(lambda kmer: len(kmers['inner_kmers'][kmer]['loci']) > 1, kmers['inner_kmers']))
        n = 0
        for kmer in kmers['junction_kmers']:
            if kmers['junction_kmers'][kmer]['source'] == 'assembly':
                continue
            if is_kmer_low_entropy(kmer):
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['junction_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['junction_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['junction_kmers'][kmer]['loci'])
            n += 1
        for kmer in kmers['inner_kmers']:
            if n > 0:
                continue
            if all_non_unique:
                continue
            if len(kmers['inner_kmers'][kmer]['loci']) > 3:
                continue
            if is_kmer_low_entropy(kmer):
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['inner_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['inner_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['inner_kmers'][kmer]['loci'])
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4)
        return path

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
                kmer['svtype'] = c.tracks[track].svtype
            kmer['coverage'] = kmer['gc_coverage']
            kmer['residue'] = 0 if 'inverse' in kmer else kmer['reference'] - r
            kmer['weight'] = 1.0
            kmer['lp_count'] = min(kmer['count'], kmer['coverage'] * kmer['reference'])
        if c.select:
            self.calculate_error_weights()

    def calculate_error_weights(self):
        for track in self.tracks:
            self.add_weights_to_track(track)

    def add_weights_to_track(self, track):
        l = self.tracks[track]['kmers']
        inner = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'inner']
        junction = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'junction']
        for kmers in [inner, junction]:
            for kmer in kmers:
                kmer['weight'] = float(len(inner) + len(junction)) / (2.0 * len(kmers))

    def generate_mps_linear_program(self):
        c = config.Configuration()
        problem = LpProblem("Nebula", LpMinimize)
        i = 0
        names = [''] * len(self.tracks)
        variables = [None] * (len(self.tracks) + 2 * len(self.lp_kmers))
        regex = re.compile('[^a-zA-Z0-9]')
        for track in self.tracks:
            name = regex.sub('_', track)
            variables[self.tracks[track]['index']] = LpVariable('c' + name, 0, 1)
            problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintLE, 'c' + name + '_ub', 1.0)
            problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintGE, 'c' + name + '_lb', 0.0)
            i += 1
        # error variables
        for index, kmer in enumerate(self.lp_kmers):
            variables[i] = LpVariable('e' + str(index))
            i += 1
        # absolute value of the error variables
        for index, kmer in enumerate(self.lp_kmers):
            variables[i] = LpVariable('l' + str(index))
            i += 1
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], kmer['weight']) for index, kmer in enumerate(self.lp_kmers)])
        problem += expr
        for i, kmer in enumerate(self.lp_kmers):
            self.add_mps_error_absolute_value_constraints(problem, variables, i)
            indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
            indices.append(len(self.tracks) + i)
            if kmer['svtype'] == 'INS' or kmer['type'] == 'junction':
                coeffs = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['lp_count'] - kmer['coverage'] * kmer['residue']
            if kmer['svtype'] == 'DEL' and kmer['type'] == 'inner':
                coeffs = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['lp_count'] - kmer['coverage'] * kmer['residue'] - sum(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
            expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
            problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
        return problem, variables

    def round_genotype(self, c, svtype):
        if c > 0.75:
            return (1.0, '1/1')
        elif c > 0.25 or c > 0.15 and svtype == 'INS':
            return (0.5, '1/0')
        else:
            return (0.0, '0/0')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcInnerKmersIntegerProgrammingJob(CgcIntegerProgrammingJob):

    _name = 'CgcInnerKmersIntegerProgrammingJob'
    _category = 'genotyping'
    _previous_job = CgcCounterJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcInnerKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    def transform(self, kmers, track_name):
        c = config.Configuration()
        lp_kmers = {}
        all_non_unique = all(map(lambda kmer: len(kmers['inner_kmers'][kmer]['loci']) > 1, kmers['inner_kmers']))
        if kmers['junction_kmers']:
            return None
        for kmer in kmers['inner_kmers']:
            if all_non_unique:
                continue
            if len(kmers['inner_kmers'][kmer]['loci']) > 3:
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['inner_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['inner_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['inner_kmers'][kmer]['loci'])
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4, sort_keys = True)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcJunctionKmersIntegerProgrammingJob(CgcIntegerProgrammingJob):

    _name = 'CgcJunctionKmersIntegerProgrammingJob'
    _category = 'preprocessing'
    _previous_job = CgcCounterJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcJunctionKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    def transform(self, kmers, track_name):
        c = config.Configuration()
        lp_kmers = {}
        all_non_unique = all(map(lambda kmer: len(kmers['inner_kmers'][kmer]['loci']) > 1, kmers['inner_kmers']))
        for kmer in kmers['junction_kmers']:
            if kmers['junction_kmers'][kmer]['source'] == 'assembly':
                continue
            if len(kmers['junction_kmers'][kmer]['loci']) > 1:
                continue
            if is_kmer_low_entropy(kmer):
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['junction_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['junction_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['junction_kmers'][kmer]['loci'])
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4, sort_keys = True)
        return path

    def round_genotype(self, c, svtype):
        if c > 0.75:
            return (1.0, '1/1')
        elif c > 0.25:#or svtype == 'INS' and c > 0.15:
            return (0.5, '1/0')
        else:
            return (0.0, '0/0')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExportGenotypingKmersJob(map_reduce.Job):

    _name = 'ExportGenotypingKmersJob'
    _category = 'genotyping'
    _previous_job = CgcIntegerProgrammingJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = ExportGenotypingKmersJob(**kwargs)
        job.execute()

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = {}
        self.load_inner_kmer_tracks()
        self.load_junction_kmer_tracks()
        self.export_kmers()
        self.export_tracks()
        exit()

    def load_junction_kmer_tracks(self):
        print('selecting junction kmers...')
        self.kmers['junction_kmers'] = {}
        job = CgcJunctionKmersIntegerProgrammingJob()
        tracks = bed.load_tracks_from_file(os.path.join(job.get_current_job_directory(), 'merge.bed'))
        self.junction_tracks = set()
        for track in tracks:
            if float(track.lp_value) > 0.25 or track.svtype == 'INS' and float(track.lp_value) > 0.15:
                kmers = json.load(open(os.path.join(job.get_current_job_directory(), track.id + '.json'), 'r'))
                self.tracks[track.id] = track
                self.junction_tracks.add(track.id)
                for kmer in kmers:
                    self.kmers['junction_kmers'][kmer] = kmers[kmer]
                    self.kmers['junction_kmers'][kmer]['count'] = 0
                    self.kmers['junction_kmers'][kmer]['total'] = 0
                    self.kmers['junction_kmers'][kmer].pop('lp_count')

    def load_inner_kmer_tracks(self):
        print('selecting inner kmers...')
        self.kmers['inner_kmers'] = {}
        job = CgcInnerKmersIntegerProgrammingJob()
        tracks = bed.load_tracks_from_file(os.path.join(job.get_current_job_directory(), 'merge.bed'))
        self.inner_tracks = set()
        for track in tracks:
            if float(track.lp_value) > 0.25:
                kmers = json.load(open(os.path.join(job.get_current_job_directory(), track.id + '.json'), 'r'))
                self.tracks[track.id] = track
                self.inner_tracks.add(track.id)
                for kmer in kmers:
                    self.kmers['inner_kmers'][kmer] = kmers[kmer]
                    self.kmers['inner_kmers'][kmer]['count'] = 0
                    self.kmers['inner_kmers'][kmer]['total'] = 0
                    self.kmers['inner_kmers'][kmer].pop('lp_count')

    def export_kmers(self):
        self._previous_job = CgcCounterJob
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
        self.kmers['gc_kmers'] = kmers['gc_kmers']
        self.kmers['depth_kmers'] = kmers['depth_kmers']
        for kmer in self.kmers['gc_kmers']:
            self.kmers['gc_kmers'][kmer]['count'] = 0
        for kmer in self.kmers['depth_kmers']:
            self.kmers['depth_kmers'][kmer]['count'] = 0
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)

    def export_tracks(self):
        with open(os.path.join(self.get_current_job_directory(), 'all.bed'), 'w') as bed_file:
            for track in self.tracks:
                bed_file.write(self.tracks[track].serialize())
        with open(os.path.join(self.get_current_job_directory(), 'both.bed'), 'w') as bed_file:
            for track in self.inner_tracks.intersection(self.junction_tracks):
                bed_file.write(self.tracks[track].serialize())
        with open(os.path.join(self.get_current_job_directory(), 'inner.bed'), 'w') as bed_file:
            for track in self.inner_tracks.difference(self.junction_tracks):
                bed_file.write(self.tracks[track].serialize())
        with open(os.path.join(self.get_current_job_directory(), 'junction.bed'), 'w') as bed_file:
            for track in self.junction_tracks.difference(self.inner_tracks):
                bed_file.write(self.tracks[track].serialize())

