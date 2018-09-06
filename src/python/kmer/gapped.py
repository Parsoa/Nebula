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
        job = CountUniqueGappedKmersJob(job_name = 'CountUniqueGappedKmersJob_', previous_job_name = 'UniqueGappedKmersJob_', category = 'programming', **kwargs)
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
                with open(os.path.join(self.get_current_job_directory(), 'unique_gapped_kmers_' + track + '.json'), 'w') as track_file:
                    json.dump(kmers, track_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            t = {}
            for track in tracks:
                t[track] = os.path.join(self.get_current_job_directory(), 'unique_gapped_kmers_' + track + '.json')
            json.dump(t, json_file, indent = 4, sort_keys = True)
        self.round_robin()

    def transform(self):
        c = config.Configuration()
        for read, name in self.parse_fastq():
            kmers = extract_canonical_gapped_kmers(read)
            for kmer in kmers:
                if kmer in self.kmers:
                    self.kmers[kmer]['count'] += 1

    def reduce(self):
        kmers = self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

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
        job = GappedKmersIntegerProgrammingJob(job_name = 'GappedKmersIntegerProgrammingJob_', previous_job_name = 'CountUniqueGappedKmersJob_', category = 'programming', batch_file_prefix = 'gapped_kmers', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        tracks = self.load_previous_job_results()
        self.round_robin(tracks)
        print(self.get_previous_job_directory())
        if c.simulation:
            self.counts_provider = counttable.DictionaryCountsProvider(json.load(open(os.path.join(self.get_previous_job_directory(), 'kmers.json'))))
        else:
            self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        self.load_reference_counts_provider()
        self.gapped_kmers = {}

    def transform(self, track, track_name):
        with open(track, 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers['inner'].keys():
                if kmer in kmers['outer']:
                    kmers['inner'].pop(kmer, None)
                    kmers['outer'].pop(kmer, None)
            for side in ['inner', 'outer']:
                for kmer in kmers[side]:
                    if not kmer in self.gapped_kmers:
                        cc = self.counts_provider.get_kmer_count(str(kmer))
                        if cc:
                            self.gapped_kmers[kmer] = {
                                'side': side,
                                'type': 'gapped',
                                'count': cc,
                                'tracks': {},
                                'reference': 1,
                            }
                    if kmer in self.gapped_kmers:
                        self.gapped_kmers[kmer]['tracks'][track_name] = kmers[side][kmer]['track']
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner': {kmer: self.gapped_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.gapped_kmers, kmers['inner']))},
                    'outer': {kmer: self.gapped_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.gapped_kmers, kmers['outer']))},
                }, json_file, indent = 4, sort_keys = True)
        return path

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + str(self.index) + '.json'), 'w')
        json.dump(self.gapped_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        c = config.Configuration()
        self.index_kmers()
        self.index_tracks()
        self.calculate_residual_coverage()
        print('exporting kmers...')
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        print('generating linear program...')
        self.solve()

    def index_kmers(self):
        c = config.Configuration()
        self.tracks = {}
        self.kmers = []
        index = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + str(i) + '.json')
            if not os.path.isfile(path):
                continue
            print(path)
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers:
                    if not kmer in index:
                        index[kmer] = len(self.kmers)
                        self.kmers.append(copy.deepcopy(kmers[kmer]))
                        self.kmers[len(self.kmers) - 1]['kmer'] = kmer
                        for track in kmers[kmer]['tracks']:
                            self.tracks[track] = True
        print(green(len(self.kmers)), 'kmers')
        return self.kmers

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.kmers:
            kmer['residue'] = 0
            kmer['coverage'] = 42#c.coverage

    def generate_linear_program(self):
        c = config.Configuration()
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        self.incorporate_gapped_kmers(problem)
        return problem

    def incorporate_gapped_kmers(self, problem):
        # the coverage of each event
        for track in self.tracks:
            tokens = track.split('_')
            problem.variables.add(names = ['c' + str(tokens[1])],
                ub = [1.0],
            )
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.kmers)],
            lb = [(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in self.kmers]
        )
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.kmers)],
            obj = [1.0] * len(self.kmers),
        )
        # constraints
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.kmers):
            if kmer['side'] == 'outer':
                # (1 - T)xR + E = C -> -TxR + E = C - R
                ind = list(map(lambda track: self.tracks[track], kmer['tracks']))
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
                ind = list(map(lambda track: self.tracks[track], kmer['tracks']))
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
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = [len(self.tracks) + len(self.kmers) + index, len(self.tracks) + index],
                    val = [1.0, 1.0],
                )],
                rhs = [0],
                senses = ['G']
            )
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = [len(self.tracks) + len(self.kmers) + index, len(self.tracks) + index],
                    val = [1.0, -1.0],
                )],
                rhs = [0],
                senses = ['G']
            )
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(self.kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
        return problem

    def solve(self):
        problem = self.generate_linear_program()
        problem.write(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_program.lp'))
        problem.solve()
        solution = problem.solution.get_values()
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_solution.json'), 'w') as json_file:
            json.dump({'variables': problem.solution.get_values()}, json_file, indent = 4, sort_keys = True)
        obj = 0
        #for i in range(len(self.tracks), len(self.tracks) + len(self.inner_kmers)):
        #    obj += abs(solution[i])
        #max_error = sum(list(map(lambda kmer: max(abs(kmer['count'] - kmer['coverage'] * kmer['residue']), abs(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks']))), self.kmers)))
        #print('error ratio:', float(obj) / max_error)
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_merge.bed'), 'w') as bed_file:
            for track in self.tracks:
                tokens = track.split('_')
                index = self.tracks[track]
                s = int(round(2 * solution[index]))
                s = '(0, 0)' if s == 2 else '(1, 0)' if s == 1 else '(1, 1)'
                bed_file.write(tokens[0] + '\t' + #0
                            tokens[1] + '\t' + #1
                            tokens[2] + '\t' + #2
                            s + '\t' + #3
                            str(solution[index]) + '\t' + #4
                            #str(len(track['inner_kmers'])) + '\t' + #5
                            self.batch_file_prefix + '\n') #6

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
