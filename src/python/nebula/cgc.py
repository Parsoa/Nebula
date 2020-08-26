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
import functools
import traceback

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
        print('Loading kmers from', green(path))
        with open(path, 'r') as json_file:
            self.kmers = json.load(json_file)
        self.gc_kmers = self.kmers['gc_kmers']
        self.depth_kmers = self.kmers['depth_kmers']
        self.inner_kmers = self.kmers['inner_kmers']
        self.junction_kmers = self.kmers['junction_kmers']
        for path in c.kmers[1:]:
            print('Loading kmers from', green(path))
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
                                    x[kmer]['loci'][locus] = kmers[kmer_type][kmer]['loci'][locus]
                                else:
                                    x[kmer]['loci'][locus]['masks'].update(kmers[kmer_type][kmer]['loci'][locus]['masks'])
                            x[kmer]['tracks'].update(kmers[kmer_type][kmer]['tracks'])
        l = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        for kmer in self.junction_kmers:
            i = len(self.junction_kmers[kmer]['loci'])
            if i < len(l):
                l[i] += 1
        for i, j in enumerate(l):
            print(i, j)
        user_print('GC kmers:', blue(len(self.gc_kmers)))
        user_print('Depth kmers:', blue(len(self.depth_kmers)))
        user_print('Inner kmers:', blue(len(self.inner_kmers)))
        user_print('Junction kmers:', blue(len(self.junction_kmers)))
        user_print('Counting..')

    def export_counter_input(self):
        for kmer in self.gc_kmers:
            self.gc_kmers[kmer]['count'] = 0
        _kmers = {}
        for kmers in [self.inner_kmers, self.junction_kmers]:
            for kmer in kmers:
                kmers[kmer]['total'] = 0
                _kmers[kmer] = {}
                _kmers[kmer]['loci'] = {}
                _kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                for locus in kmers[kmer]['loci']:
                    _kmers[kmer]['loci'][locus] = {}
                    _kmers[kmer]['loci'][locus]['masks'] = kmers[kmer]['loci'][locus]['masks']
        for kmer in self.gc_kmers:
            _kmers[kmer] = {}
            _kmers[kmer]['loci'] = {}
        for kmer in self.depth_kmers:
            _kmers[kmer] = {}
            _kmers[kmer]['loci'] = {}
        if self.resume_from_reduce:
            return
        system_print('Dumping', green(len(_kmers)), 'kmers for the counter..')
        with open(os.path.join(self.get_current_job_directory(), 'raw_kmers.json'), 'w') as json_file:
            json.dump(_kmers, json_file, indent = 4)

    def merge_count(self, kmer, tokens):
        count = tokens[0] 
        total = tokens[1] 
        canon = canonicalize(kmer)
        if canon in self.gc_kmers:
            self.gc_kmers[canon]['count'] += count
        if canon in self.depth_kmers:
            self.depth_kmers[canon]['count'] += count
        if canon in self.inner_kmers:
            self.inner_kmers[canon]['count'] += count
            self.inner_kmers[canon]['total'] += total
        if canon in self.junction_kmers:
            self.junction_kmers[canon]['count'] += count
            self.junction_kmers[canon]['total'] += total

    def reduce(self):
        self.merge_counts()
        self.export_counted_kmers()
        self.stats = self.estimate_depth_of_coverage()
        self.estimate_gc_content_coverage()
        self.export_genotyping_tracks()
        return self.stats

    def estimate_depth_of_coverage(self):
        c = config.Configuration()
        counts = list(map(lambda kmer: self.depth_kmers[kmer]['count'], self.depth_kmers))
        print('Estimating coverage statisics from', len(counts), 'depth kmers..')
        visualizer.histogram(counts, 'depth_unfiltered', self.get_current_job_directory(), 'count', 'kmers', 1)
        while True:
            std = np.std(counts)
            mean = np.mean(counts)
            counts = list(filter(lambda x: x <= 3 * mean, counts))
            if abs(np.mean(counts) - mean) < 0.5:
                break
        print('Distribution:', str(int(mean)) + 'x with STD of', str(int(std)) + '.')
        #
        visualizer.histogram(counts, 'depth_filtered', self.get_current_job_directory(), 'count', 'kmers', 1)
        stats = {'coverage': int(round(mean)), 'std': std}
        with open(os.path.join(self.get_current_job_directory(), 'stats_' + str(c.ksize) + '.json'), 'w') as json_file:
            json.dump(stats, json_file, sort_keys = True, indent = 4)
        return stats

    def estimate_gc_content_coverage(self):
        c = config.Configuration()
        print('Estimating coverage per GC content level..')
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

    def adjust_kmer_gc_coverage(self, kmer):
        kmer['gc_coverage'] = sum(map(lambda locus: self.gc[int(kmer['loci'][locus]['gc'])], kmer['loci'])) // len(kmer['loci'])

    def export_counted_kmers(self):
        c = config.Configuration()
        print('Exporting kmers..')
        if self.resume_from_reduce:
            return
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
            self.adjust_kmer_gc_coverage(self.inner_kmers[kmer])
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.junction_kmers:
            self.junction_kmers[kmer]['gc_coverage'] = self.stats['coverage']
            if self.junction_kmers[kmer]['inverse']:
                self.junction_kmers[kmer]['count'] = self.junction_kmers[kmer]['total'] - self.junction_kmers[kmer]['count']
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        print('Exporting kmers for', len(self.tracks), 'tracks.')
        exporter = map_reduce.TrackExportHelper(tracks = self.tracks, output_dir = self.get_current_job_directory(), num_threads = 12)
        exporter.execute()
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: os.path.join(self.get_current_job_directory(), track + '.json') for track in self.tracks}, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'all_tracks.bed'), 'w') as bed_file:
            bed_file.write('\t'.join([str(s) for s in ['#CHROM', 'BEGIN', 'END', 'KMERS', 'GENOTYPE']]) + '\n')
            for track in self.tracks:
                if track in c.tracks:
                    t = c.tracks[track]
                    l = len(self.tracks[track]['inner_kmers']) + len(self.tracks[track]['junction_kmers'])
                    bed_file.write('\t'.join([str(s) for s in [t.chrom, t.begin, t.end, l, t.genotype]]) + '\n')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class PreFilteringJob(map_reduce.Job):

    _name = 'PreFilteringJob'

    def load_inputs(self):
        c = config.Configuration()
        self.tracks = []
        for i, bed_path in enumerate(c.bed):
            self.tracks.append(bed.as_dict(bed.sort_tracks(bed.filter_short_tracks(bed.load_tracks_from_file(bed_path)))))
        self.all_tracks = json.load(open(os.path.join(c.workdir, c.samples[0], 'CgcCounterJob', 'batch_merge.json')))
        self.round_robin(self.all_tracks)

    def calc_kmer_genotype_likelihood(self, kmer, genotype, std):
        c = 0 if genotype == '0/0' else 0.5 if genotype == '1/0' else 1.0
        s = 0.25 if genotype == '0/0' else 0.50 if genotype == '1/0' else 1.00
        if kmer['type'] == 'junction':
            d = statistics.NormalDistribution(kmer['gc_coverage'] * c, std * s)
        else:
            if kmer['trend'] == 'upward':
                d = statistics.NormalDistribution(kmer['gc_coverage'] * c, std * s)
            else:
                d = statistics.NormalDistribution(kmer['gc_coverage'] * (1 - c), std * (0.25 / s))
        return d.log_pmf(kmer['count'])

    def transform(self, _, track):
        if not '39061957' in track:
            return None
        print(track)
        c = config.Configuration()
        genotypes = ['0/0', '1/0', '1/1']
        kmers = {}
        scores = {}
        for i, sample in enumerate(c.samples):
            kmers[sample] = json.load(open(os.path.join(c.workdir, sample, 'CgcCounterJob', track + '.json')))
            for kmer_type in ['junction_kmers', 'inner_kmers']:
                for kmer in kmers[sample][kmer_type]:
                    #if kmer != 'AAAAAGAAAGAAAGAAAGAAAGATAAATGTCT':
                    if kmer != 'CCCCGCAACCCCAACCGGGGCGCAAGAGAAAC':
                        continue
                    if len(kmers[sample][kmer_type][kmer]['loci']) > 1:
                        continue
                    if len(kmers[sample][kmer_type][kmer]['tracks']) > 1:
                        continue
                    if i == 0:
                        scores[kmer] = []
                    likelihoods = {g: self.calc_kmer_genotype_likelihood(kmers[sample][kmer_type][kmer], g, c.std) for g in genotypes}
                    a = self.tracks[i][track].genotype if track in self.tracks[i] else '0/0'
                    other_likelihoods = {g: likelihoods[g] for g in genotypes if g != a}
                    g = max(other_likelihoods, key = other_likelihoods.get)
                    s = likelihoods[a] - likelihoods[g]
                    print(sample)
                    print(kmers[sample][kmer_type][kmer]['count'], likelihoods, g, a, s)
                    scores[kmer].append(s)
        remove = {}
        for kmer in scores:
            #score = functools.reduce(lambda x, y: x * y, scores[kmer])
            score = sum(scores[kmer])
            #print(score)
            if score >= 0.0:
                #print(cyan('Keeping', kmer, scores[kmer]))
                pass
            else:
                #print(yellow('Filtering', kmer, scores[kmer]))
                remove[kmer] = True
        debug_breakpoint()
        for kmer in remove:
            scores.pop(kmer)
        if len(scores) > 0:
            for i, sample in enumerate(c.samples):
                path = os.path.join(c.workdir, sample, 'PreFilteringJob')
                if not os.path.exists(path):
                    os.makedirs(path)
                with open(os.path.join(path, track + '.json'), 'w') as json_file:
                    _kmers = {}
                    for kmer_type in ['junction_kmers', 'inner_kmers']:
                        _kmers[kmer_type] = {}
                        for kmer in kmers[sample][kmer_type]:
                            if kmer in scores:
                                _kmers[kmer_type][kmer] = kmers[sample][kmer_type][kmer]
                    json.dump(_kmers, json_file, indent = 4)
            return track
        else:
            return None
    
    def reduce(self):
        c = config.Configuration()
        for sample in c.samples:
            path = os.path.join(c.workdir, sample, 'PreFilteringJob')
            tracks = {}
            for batch in self.load_output():
                for track in batch:
                    tracks[track] = os.path.join(path, track + '.json')
            with open(os.path.join(path, 'batch_merge.json'), 'w') as json_file:
                json.dump(tracks, json_file, indent = 4)

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
    def transform(self, path, track_name):
        c = config.Configuration()
        kmers = json.load(open(path))
        lp_kmers = {}
        genotypes = ['0/0', '1/0', '1/1']
        all_non_unique = all(map(lambda kmer: len(kmers['inner_kmers'][kmer]['loci']) > 1, kmers['inner_kmers']))
        n = 0
        for kmer_type in ['junction_kmers', 'inner_kmers']:
            for kmer in kmers[kmer_type]:
                if len(kmers[kmer_type][kmer]['loci']) > 1:
                    continue
                if len(kmers[kmer_type][kmer]['tracks']) > 1:
                    continue
                self.calculate_residual_coverage(kmers[kmer_type][kmer])
                if c.select:
                    p = False
                    for track in kmers[kmer_type][kmer]['tracks']:
                        t = bed.track_from_id(track_name)
                        likelihoods = {g: self.calc_kmer_genotype_likelihood(kmers[kmer_type][kmer], g, c.std) for g in genotypes}
                        g = max(likelihoods, key = likelihoods.get)
                        if '1' in g:
                            if track in c.tracks and '1' in c.tracks[track].genotype:
                                p = True
                        else:
                            if not track in c.tracks or not '1' in c.tracks[track].genotype:
                                p = True
                    if not p:
                        continue
                lp_kmers[kmer] = True
                self.lp_kmers[kmer] = kmers[kmer_type][kmer]
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        return path

    def calculate_residual_coverage(self, kmer):
        r = 0
        kmer['reduction'] = len(kmer['loci'])
        for track in kmer['tracks']:
            r += kmer['tracks'][track]
        kmer['weight'] = 1.0
        kmer['residue'] = 0 if 'inverse' in kmer else kmer['reduction'] - r
        kmer['coverage'] = kmer['gc_coverage']
        kmer['lp_count'] = min(kmer['count'], kmer['coverage'] * kmer['reduction'])

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
            if kmer['type'] == 'junction' or kmer['trend'] == 'upward':
                coeffs = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['lp_count'] - kmer['coverage'] * kmer['residue']
            else:
                coeffs = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['lp_count'] - kmer['coverage'] * kmer['residue'] - sum(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
            expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
            problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
        return problem, variables

    def round_genotype(self, c, svtype):
        if c >= 0.75:
            return (1.0, '1/1')
        elif c >= 0.25:
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
        self.load_kmers()
        self.export_kmers()
        self.export_tracks()
        exit()

    def load_kmers(self):
        c = config.Configuration()
        print('loading kmers..')
        self.kmers['inner_kmers'] = {}
        self.kmers['junction_kmers'] = {}
        job = CgcIntegerProgrammingJob()
        tracks = bed.load_tracks_from_file(os.path.join(job.get_current_job_directory(), 'merge.bed'))
        print('Loaded', len(tracks), 'tracks.')
        for track in tracks:
            if (track.genotype == '0/0' and track.id not in c.tracks) or ((track.id in c.tracks) and (('1' in track.genotype and '1' in c.tracks[track.id].genotype) or ('1' not in track.genotype and '1' not in c.tracks[track.id].genotype))):
                self.tracks[track.id] = track
                kmers = json.load(open(os.path.join(job.get_current_job_directory(), track.id + '.json'), 'r'))
                for kmer in kmers:
                    kmer_type = kmers[kmer]['type'] + '_kmers'
                    self.kmers[kmer_type][kmer] = kmers[kmer]
                    self.kmers[kmer_type][kmer]['count'] = 0
                    self.kmers[kmer_type][kmer]['total'] = 0
                    self.kmers[kmer_type][kmer].pop('lp_count')

    def export_kmers(self):
        self._previous_job = CgcCounterJob
        print('Exporting kmers for genotyping:')
        print(len(self.kmers['inner_kmers']), 'inner kmers.')
        print(len(self.kmers['junction_kmers']), 'junction kmers.')
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

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class MultiSampleExportGenotypingKmersJob(map_reduce.Job):

    _name = 'ExportGenotypingKmersJob'
    _category = 'genotyping'
    _previous_job = CgcIntegerProgrammingJob

    def load_inputs(self):
        c = config.Configuration()
        kmers = {}
        calls = []
        tracks = []
        inner_kmers = {}
        junction_kmers = {}
        for i, bed_path in enumerate(c.bed):
            path = os.path.join(self.get_output_directory(), c.samples[i], 'CgcIntegerProgrammingJob', 'merge.bed')
            calls.append(bed.load_tracks_from_file_as_dict(path))
            tracks.append(bed.as_dict(bed.sort_tracks(bed.filter_short_tracks(bed.load_tracks_from_file(bed_path)))))
        m = 0
        k = 0
        u = 0
        print(len(c.tracks), 'tracks in total.')
        for track in c.tracks:
            # one sample has this event as present
            if any(['1' in tracks[i][track].genotype for i in range(len(tracks)) if track in tracks[i]]):
                m += 1
                if any(['1' in calls[i][track].genotype and '1' in tracks[i][track].genotype for i in range(len(tracks)) if track in calls[i] and track in tracks[i]]):
                    k += 1
                    found = False
                    for i, sample in enumerate(c.samples):
                        if track in calls[i]:
                            if not found:
                                u += 1
                                found = True
                            g_1 = calls[i][track].genotype
                            if (track in tracks[i] and (g_1 == tracks[i][track].genotype or '1' in g_1 and '1' in tracks[i][track].genotype)) or (track not in tracks[i] and not '1' in g_1):
                            #if (track in tracks[i] and (g_1 == tracks[i][track].genotype)) or (track not in tracks[i] and not '1' in g_1):
                                with open(os.path.join(self.get_output_directory(), sample, 'CgcIntegerProgrammingJob', str(track) + '.json'), 'r') as json_file:
                                    _kmers = json.load(json_file)
                                if not track in kmers:
                                    kmers[track] = _kmers
                                else:
                                    kmers[track] = {kmer: kmers[track][kmer] for kmer in kmers[track] if kmer in _kmers}
                    if track in kmers:
                        inner_kmers.update({kmer: kmers[track][kmer] for kmer in kmers[track] if kmers[track][kmer]['type'] == 'inner'})
                        junction_kmers.update({kmer: kmers[track][kmer] for kmer in kmers[track] if kmers[track][kmer]['type'] == 'junction'})
        print(m)
        print(u)
        print(k)
        n = len([track for track in kmers if track in kmers and len(kmers[track]) != 0])
        print(n, 'tracks.')
        print(len(inner_kmers), 'inner kmers.')
        print(len(junction_kmers), 'junction kmers.')
        with open(os.path.join(self.get_output_directory(), sample, 'CgcCounterJob', 'kmers.json'), 'r') as json_file:
            _kmers = json.load(json_file)
            gc_kmers = _kmers['gc_kmers']
            depth_kmers = _kmers['depth_kmers']
        for _kmers in [junction_kmers, inner_kmers, depth_kmers, gc_kmers]:
            for kmer in _kmers:
                _kmers[kmer]['count'] = 0
                _kmers[kmer]['total'] = 0
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump({'junction_kmers': junction_kmers, 'inner_kmers': inner_kmers, 'depth_kmers': depth_kmers, 'gc_kmers': gc_kmers}, json_file, indent = 4)
        exit()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class PcaClusteringJob(map_reduce.Job):

    _name = 'PcaClusteringJob'
    _category = 'misc'
    _previous_job = None 

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    def load_inputs(self):
        from sklearn.decomposition import PCA
        from matplotlib import pyplot as plt
        c = config.Configuration()
        payload = json.load(open(os.path.join(c.workdir, 'simons.json'), 'r'))
        regions = json.load(open(os.path.join(c.workdir, 'region.json'), 'r'))
        features = []
        populations = []
        for name in payload['data']:
            if name in regions:
                #if regions[name] == 'Africa':
                #    continue
                self.tracks = bed.load_tracks_from_file(os.path.join(c.workdir, 'Genotypes', name + '.bed' ))
                features.append([self.ordinalize_genotype(t.lp_genotype) for t in self.tracks if 'INS' in t.id])
                populations.append(regions[name])

        print(len(features))
        print(len(populations))

        #pca = PCA(n_components = 2)
        #components = pca.fit_transform(features)

        #print(sorted(set(populations)))

        #fig = plt.figure(figsize = (8,8))
        #ax = fig.add_subplot(1,1,1)
        #ax.set_xlabel('Principal Component 1', fontsize = 15)
        #ax.set_ylabel('Principal Component 2', fontsize = 15)
        #ax.set_title('2 component PCA', fontsize = 20)
        ##colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        #colors = ['xkcd:' + a for a in ['blue', 'green', 'red', 'pink', 'brown', 'purple', 'yellow']]
        #for j, p in enumerate(sorted(set(populations))):
        #    print(p, colors[j])
        #    a = [f[0] for i, f in enumerate(components) if populations[i] == p]
        #    b = [f[1] for i, f in enumerate(components) if populations[i] == p]
        #    print(a)
        #    print(b)
        #    print('====')
        #    ax.scatter(a, b, color = colors[j])
        #ax.legend(sorted(set(populations)))
        ##fig.savefig(os.path.join(c.workdir, 'pca.png'))
        #plt.close(fig)
        #  Hardy-Weinberg equlibrium
        heteros = [0 for f in features[0]]
        genotypes = [[0, 0, 0] for f in features[0]]
        genotype_counts = [[0, 0, 0] for f in features[0]]
        alleles = [0 for f in features[0]]
        for feature in features:
            for i, f in enumerate(feature):
                alleles[i] += f #1 if f != 0 else 0
                heteros[i] += 1 if f == 1 else 0
                genotypes[i][f] += 1
                genotype_counts[i][f] += 1
        m = 0
        for i in range(len(alleles)):
            l = len(features)
            alleles[i] = alleles[i] / float(2 * len(features))
            heteros[i] = heteros[i] / float(len(features))
            if heteros[i] > 0.90 and alleles[i] > 0.5:
                print(self.tracks[i])
            genotypes[i] = [str(g / float(l)) for g in genotypes[i]]
            genotype_counts[i] = [str(g) for g in genotype_counts[i]]
        with open(os.path.join(self.get_current_job_directory(), 'x_INS.txt'), 'w') as txt_file:
            for genotype in genotypes:
                txt_file.write(' '.join(genotype) + '\n')
        with open(os.path.join(self.get_current_job_directory(), 'c_INS.txt'), 'w') as txt_file:
            for genotype_count in genotype_counts:
                txt_file.write(' '.join(genotype_count) + '\n')
        visualizer.scatter(alleles, heteros, 'Hardy Weinberg principle', self.get_current_job_directory(), 'Allele Frequency', 'Genotype Frequency')
        exit()

    def ordinalize_genotype(self, g):
        if g == '0/0':
            return 0
        if g == '1/0' or g == '0/1':
            return 1
        return 2
