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

from kmer import (
    bed,
    sets,
    config,
    counter,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
)

from kmer.sv import StructuralVariation, Inversion, Deletion, SNP
from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import cplex
import numpy
import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job for finding StructuralVariation breakpoints
# Algorithm: starts with a set of structural variation events and their approximate breakpoints and tries to refine them
# considers a radius of [-50, 50] around each end of the breakpoint, and for each pair of endpoints within that radius considers
# the area as the structural variations and applies it to the reference genome to generate a set of kmers. Discards those endpoints
# whose kmers do not all appear in the base genome the event was detected in. 
# Output: Reports the remaining boundary candidates with their list of associated kmers and the count of those kmers in the
# base genome.
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractGappedKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractGappedKmersJob(job_name = 'ExtractGappedKmersJob_', previous_job_name = '', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.bedtools = {str(track): track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'all.bed'))}
        self.round_robin(self.bedtools, lambda track: re.sub(r'\s+', '_', str(track).strip()).strip(), lambda track: track.end - track.start > 1000000) 

    def transform(self, track, track_name):
        sv = self.get_sv_type()(track)
        c = config.Configuration()
        gapped_kmers = sv.extract_boundary_gapped_kmers()
        path = os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump({'gapped_kmers': gapped_kmers}, json_file, sort_keys = True, indent = 4)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueGappedKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueGappedKmersJob(job_name = 'UniqueGappedKmersJob_', previous_job_name = 'ExtractGappedKmersJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        #self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.gapped_counts_provider = counttable.JellyfishCountsProvider(self.get_gapped_reference_counts_provider())
        self.gapped_kmers = {'inner': {}, 'outer': {}}
        tracks = self.load_previous_job_results()
        for track in tracks:
            with open(tracks[track]) as json_file:
                gapped_kmers = json.load(json_file)['gapped_kmers']
                for side in self.gapped_kmers:
                    for kmer in gapped_kmers[side]:
                        k = canonicalize(kmer)
                        k = k[:15] + k[20:]
                        if not k in self.gapped_kmers[side]:
                            self.gapped_kmers[side][k] = {}
                        if track in self.gapped_kmers[side][k]:
                            self.gapped_kmers[side][k][track] += 1
                        else:
                            self.gapped_kmers[side][k][track] = 1
        self.round_robin(tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(track, 'r') as json_file:
            _gapped_kmers = json.load(json_file)['gapped_kmers']
        gapped_kmers = {'inner': {}, 'outer': {}}
        n = 0
        for kmer in _gapped_kmers['inner'].keys():
            if kmer in _gapped_kmers['outer']:
                _gapped_kmers['inner'].pop(kmer, None)
                _gapped_kmers['outer'].pop(kmer, None)
        for side in ['inner', 'outer']:
            for kmer in _gapped_kmers[side]:
                k = canonicalize(kmer)
                k = k[:15] + k[20:]
                #if kmer != self.gapped_kmers[side][k][track_name]:
                #    print(yellow(kmer), self.gapped_kmers[side][k][track_name])
                s = self.get_novelty_score(kmer)
                #count = self.gapped_counts_provider.get_kmer_count(kmer)
                gapped_kmers[side][k] = {'unique': len(self.gapped_kmers[side][k]) == 1, 'score': s, 'track': self.gapped_kmers[side][k][track_name] }
                n += 1
                print(n, 'of', len(gapped_kmers['inner']) + len(gapped_kmers['outer']))
        path = os.path.join(self.get_current_job_directory(), 'unique_gapped_kmers_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump({'gapped_kmers': gapped_kmers}, json_file, sort_keys = True, indent = 4)
        return path

    def get_novelty_score(self, kmer):
        l = kmer[:15]
        r = kmer[20:]
        #print(kmer)
        #print(green(l) + white(kmer[15:20]) + green(r))
        counts = {}
        self.recurse(l, r, ['A', 'T', 'C', 'G'], counts, kmer)
        s = list(map(lambda x: counts[x], counts))
        return sum(s)

    def recurse(self, left, right, bases, counts, kmer):
        if len(left + right) == 35:
            #print(blue(left[:15]) + white(left[15:]) + blue(right))
            c = self.gapped_counts_provider.get_kmer_count(left + right)
            counts[left + right] = c
            return
        for base in bases:
            self.recurse(left + base, right, bases, counts, kmer)

    def plot(self, tracks):
        x = {'inner': [], 'outer': []}
        y = {'inner': [], 'outer': []}
        print(len(tracks))
        for track in tracks:
            print(track)
            with open(tracks[track], 'r') as json_file:
                gapped_kmers = json.load(json_file)['gapped_kmers']
                for side in ['inner', 'outer']:
                    print(side)
                    x[side].append(len(list(filter(lambda kmer: gapped_kmers[side][kmer]['unique'], gapped_kmers[side]))))
                    y[side] += list(map(lambda kmer: gapped_kmers[side][kmer]['score'], gapped_kmers[side]))
        print(len(x['inner']))
        print(len(x['outer']))
        print(x['inner'])
        #visualizer.histogram(x = x['inner'], name = 'num unique inner gapped kmers', x_label = 'number of unique kmers', y_label = 'number of events', path = self.get_current_job_directory())
        #visualizer.histogram(x = x['outer'], name = 'num unique outer gapped kmers', x_label = 'number of unique kmers', y_label = 'number of events', path = self.get_current_job_directory())
        visualizer.histogram(x = y['inner'], name = 'number of unique inner isomers', x_label = 'percentage of novel or unique isomers', y_label = 'number of kmers', path = self.get_current_job_directory())
        visualizer.histogram(x = y['outer'], name = 'number of novel outer isomers', x_label = 'percentage of novel or unique isomers', y_label = 'number of kmers', path = self.get_current_job_directory())

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountUniqueGappedKmersJob(counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountUniqueGappedKmersJob(job_name = 'CountUniqueGappedKmersJob_', previous_job_name = 'UniqueGappedKmersJob_', category = 'programming', kmer_type = 'unique_gapped', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        tracks = self.load_previous_job_results()
        for track in tracks:
            print(cyan(track))
            with open(tracks[track], 'r') as json_file:
                kmers = json.load(json_file)['gapped_kmers']
                for kmer in kmers['inner']:
                    # is unique to this event and none of it's isomers ever appear
                    if kmers['inner'][kmer]['unique'] and kmers['inner'][kmer]['score'] == kmers['inner'][kmer]['track']:
                        k = canonicalize(kmer)
                        self.kmers[k] = {'track': track, 'side': 'inner', 'count': 0}
                for kmer in kmers['outer']:
                    # is unique to this event and none of it's isomers appear anywhere else
                    if kmers['outer'][kmer]['unique'] and kmers['outer'][kmer]['score'] == 0:
                        k = canonicalize(kmer)
                        self.kmers[k] = {'track': track, 'side': 'outer', 'count': 0}
        self.round_robin()

    def transform(self):
        c = config.Configuration()
        for read, name in self.parse_fastq():
            kmers = extract_canonical_gapped_kmers(read)
            for kmer in kmers:
                if kmer in self.kmers:
                    self.kmers[kmer]['count'] += 1

    def reduce(self):
        self.kmers = self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        # output kmers per track
        for kmer in self.kmers:
            track = self.kmers[kmer]['track']
            if not track in self.tracks:
                self.tracks[track] = {'inner': {}, 'outer': {}}
            self.tracks[track][self.kmers[kmer]['side']][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), self.kmer_type + '_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: self.kmer_type + '_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class GappedKmersIntegerProgrammingJob(programming.IntegerProgrammingJob):

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = GappedKmersIntegerProgrammingJob(job_name = 'GappedKmersIntegerProgrammingJob_', previous_job_name = 'CountUniqueGappedKmersJob_', category = 'programming', batch_file_prefix = 'gapped_kmers', kmer_type = 'gapped', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        with open(track, 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers['inner'].keys():
                if kmer in kmers['outer']:
                    kmers['inner'].pop(kmer, None)
                    kmers['outer'].pop(kmer, None)
            for side in ['inner', 'outer']:
                for kmer in kmers[side]:
                    if not kmer in self.lp_kmers:
                        cc = self.counts_provider.get_kmer_count(str(kmer))
                        if cc:
                            self.lp_kmers[kmer] = {
                                'side': side,
                                'type': 'gapped',
                                'count': cc,
                                'tracks': {},
                                'reference': 1,
                            }
                    if kmer in self.lp_kmers:
                        self.lp_kmers[kmer]['tracks'][track_name] = kmers[side][kmer]['track']
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner': {kmer: self.lp_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.lp_kmers, kmers['inner']))},
                    'outer': {kmer: self.lp_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.lp_kmers, kmers['outer']))},
                }, json_file, indent = 4, sort_keys = True)
        return path

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            kmer['residue'] = 0
            kmer['coverage'] = 42#c.coverage

    def generate_linear_program(self):
        c = config.Configuration()
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        # the coverage of each event
        for track in self.tracks:
            tokens = track.split('_')
            problem.variables.add(names = ['c' + str(tokens[1])],
                ub = [1.0],
            )
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            lb = [(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in self.lp_kmers]
        )
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [1.0] * len(self.lp_kmers),
        )
        # constraints
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.lp_kmers):
            if kmer['side'] == 'outer':
                # (1 - T)xR + E = C -> -TxR + E = C - R
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - kmer['coverage'] * kmer['residue'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['E']
                )
            else:
                # TxR + E = C
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']],
                    senses = ['E']
                )
            self.add_error_absolute_value_constraints(problem, index)
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(self.lp_kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
        return problem

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
