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

from kmer import (
    bed,
    config,
    counter,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
)

from kmer.sv import StructuralVariation, Inversion, Deletion
from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

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

    _name = 'ExtractGappedKmersJob'
    _category = 'programming'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        extract_whole_genome()
        self.tracks = self.load_tracks()
        self.round_robin(self.tracks, filter_func = lambda track: track.end - track.begin > 1000000)

    # These kmers ARE NOT CANONICAL
    def transform(self, track, track_name):
        print(cyan(track_name))
        c = config.Configuration()
        gapped_kmers = track.extract_boundary_gapped_kmers()
        path = os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump(gapped_kmers, json_file, sort_keys = True, indent = 4)
        return 'gapped_kmers_' + track_name + '.json'

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueGappedKmersJob(map_reduce.Job):

    _name = 'UniqueGappedKmersJob'
    _category = 'programming'
    _previous_job = ExtractGappedKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        #self.gapped_counts_provider = counttable.JellyfishCountsProvider(self.get_gapped_reference_counts_provider())
        self.gapped_kmers = {'inner': {}, 'outer': {}}
        tracks = self.load_previous_job_results()
        for track in tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                gapped_kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    for kmer in gapped_kmers[side]:
                        k = canonicalize(kmer)
                        k = k[:c.hsize] + k[-c.hsize:]
                        if not k in self.gapped_kmers[side]:
                            self.gapped_kmers[side][k] = {'tracks': {}}
                        if track in self.gapped_kmers[side][k]:
                            self.gapped_kmers[side][k]['tracks'][track] += 1
                        else:
                            self.gapped_kmers[side][k]['tracks'][track] = 1
        self.round_robin(tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            _gapped_kmers = json.load(json_file)
        gapped_kmers = {'inner': {}, 'outer': {}}
        for side in ['inner', 'outer']:
            for kmer in _gapped_kmers[side]:
                if 'N' in kmer:
                    continue
                k = canonicalize(kmer)
                k = k[:c.hsize] + k[-c.hsize:]
                unique = 0
                unique += len(self.gapped_kmers['inner'][k]['tracks']) if k in self.gapped_kmers['inner'] else 0
                unique += len(self.gapped_kmers['outer'][k]['tracks']) if k in self.gapped_kmers['outer'] else 0
                # The kmer is not shared with a kmer of opposite side of another event
                if unique == 1:
                    gapped_kmers[side][k] = self.gapped_kmers[side][k]
        path = os.path.join(self.get_current_job_directory(), 'unique_gapped_kmers_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump(gapped_kmers, json_file, sort_keys = True, indent = 4)
        return 'unique_gapped_kmers_' + track_name  + '.json'

    def plot(self, tracks):
        x = {'inner': [], 'outer': []}
        y = {'inner': [], 'outer': []}
        print(len(tracks))
        for track in tracks:
            print(track)
            with open(os.path.join(self.get_current_job_directory(), tracks[track]), 'r') as json_file:
                gapped_kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    print(side)
                    x[side].append(len(gapped_kmers[side]))
        visualizer.histogram(x = x['inner'], name = 'number of unique inner gapped kmers', x_label = 'number of unique gapped kmers', y_label = 'number of events', path = self.get_current_job_directory())
        visualizer.histogram(x = x['outer'], name = 'number of unique outer gapped kmers', x_label = 'number of unique gapped kmers', y_label = 'number of events', path = self.get_current_job_directory())

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueGappedKmersScoringJob(map_reduce.Job):

    _name = 'UniqueGappedKmersScoringJob'
    _category = 'programming'
    _previous_job = UniqueGappedKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueGappedKmersScoringJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = self.load_previous_job_results()
        self.half_mers = {}
        for track in self.tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), self.tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for side in ['outer', 'inner']:
                    for kmer in kmers[side]:
                        if len(kmer) != c.ksize:
                            print(red(kmer))
                        self.kmers[kmer] = kmers[side][kmer]
                        self.kmers[kmer]['count'] = 0
                        self.kmers[kmer]['side'] = side
                        left = kmer[:c.hsize]
                        right = kmer[-c.hsize:]
                        self.half_mers[left] = True
                        self.half_mers[right] = True
                        self.half_mers[reverse_complement(left)] = True
                        self.half_mers[reverse_complement(right)] = True
        print('scoring', len(self.kmers), 'kmers')
        self.chroms = extract_whole_genome()
        for chrom in self.chroms:
            print(sys.getsizeof(self.chroms[chrom]))
        self.round_robin(self.chroms)

    def transform(self, sequence, chrom):
        c = config.Configuration()
        t = time.time()
        l = len(sequence)
        for index in range(0, l - 50):
            if index % 100000 == 1:
                s = time.time()
                p = index / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
            half_mer = sequence[index: index + c.hsize]
            if not half_mer in self.half_mers:
                continue
            for k in range(32, 43):
                kmer = canonicalize(sequence[index: index + k])
                #kmer = sequence[index: index + k]
                kmer = kmer[:c.hsize] + kmer[-c.hsize:]
                if kmer in self.kmers:
                    print(kmer)
                    if self.kmers[kmer]['side'] == 'inner' and k != c.ksize + 5:
                        continue
                    self.kmers[kmer]['count'] += 1

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        self.kmers = self.merge_counts()
        self.tracks = self.load_previous_job_results()
        for track in self.tracks.keys():
            with open(os.path.join(self.get_previous_job_directory(), self.tracks[track]), 'r') as json_file:
                self.tracks[track] = json.load(json_file)
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                self.tracks[track][self.kmers[kmer]['side']][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'gapped_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class SelectUniqueGappedKmersJob(counter.BaseExactCountingJob):

    _name = 'SelectUniqueGappedKmersJob'
    _category = 'programming'
    _previous_job = UniqueGappedKmersScoringJob
    _counter_mode = 2

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = SelectUniqueGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.half_mers = {}
        tracks = self.load_previous_job_results()
        x = []
        m = 0
        for track in tracks:
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['outer']:
                    x.append(kmers['outer'][kmer]['count'])
                    if kmers['outer'][kmer]['count'] == 0:
                        left = kmer[:c.hsize]
                        right = kmer[-c.hsize:]
                        self.kmers[kmer] = {'tracks': kmers['outer'][kmer]['tracks'], 'side': 'outer', 'count': {}}
                        for i in range(0, 11):
                            self.kmers[kmer]['count'][i] = 0
                        if not left in self.half_mers:
                            self.half_mers[left] = {}
                        self.half_mers[left][right] = kmer 
                        left = reverse_complement(left)
                        right = reverse_complement(right)
                        if not right in self.half_mers:
                            self.half_mers[right] = {}
                        self.half_mers[right][left] = kmer 
                        m += 1
        visualizer.histogram(x, 'outer_gapped_kmer_count', self.get_current_job_directory(), 'count', 'number of kmers')
        print(len(tracks), 'total events')
        print(len(self.kmers), 'kmers')
        print(m, 'events with kmers')
        self.export_accelerator_input()
        self.round_robin()

    def export_accelerator_input(self):
        with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
            json.dump(self.half_mers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)

    def reduce(self):
        c = config.Configuration()
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                line = json_file.readline()
                while line:
                    kmer = line[:line.find(':')]
                    index = line.find(':')
                    for i in range(0, 11):
                        count = int(line[index + 1: line.find(':', index + 1)])
                        self.kmers[kmer]['count'][i] += count
                        index = line.find(':', index + 1)
                    line = json_file.readline()
        self.tracks = {}
        print(len(self.kmers))
        keys = list(self.kmers.keys())
        for kmer in keys:
            side = self.kmers[kmer]['side']
            gap = max(self.kmers[kmer]['count'].items(), key = operator.itemgetter(1))[0]
            total = sum(map(lambda x: self.kmers[kmer]['count'][x], self.kmers[kmer]['count']))
            self.kmers[kmer]['gap'] = gap if self.kmers[kmer]['count'][gap] != 0 else -1
            self.kmers[kmer]['count'] = self.kmers[kmer]['count'][gap]
            self.kmers[kmer]['total'] = total
            self.kmers[kmer]['actual_gap'] = -1
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner': {}, 'outer': {}}
                self.tracks[track][side][kmer] = self.kmers[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'gapped_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'gaps.bed'), 'w') as bed_file:
            for kmer in self.kmers:
                for track in self.kmers[kmer]['tracks']:
                    t = bed.track_from_name(track)
                    bed_file.write(t.chrom + '\t' + str(t.begin) + '\t' + str(t.end) + '\t' + str(self.kmers[kmer]['gap']) + '\n')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountUniqueGappedKmersJob(map_reduce.FirstGenotypingJob, SelectUniqueGappedKmersJob):

    _name = 'CountUniqueGappedKmersJob'
    _category = 'programming'
    _previous_job = SelectUniqueGappedKmersJob
    _counter_mode = 1

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountUniqueGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        tracks = self.load_previous_job_results()
        self.half_mers = {}
        n = 0
        for track in tracks:
            n += 1
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    for kmer in kmers[side]:
                        if kmers[side][kmer]['gap'] != -1:
                            left = kmer[:c.hsize]
                            right = kmer[-c.hsize:]
                            self.kmers[kmer] = kmers[side][kmer]
                            self.kmers[kmer]['count'] = 0
                            self.kmers[kmer]['doubt'] = 0
                            print(type(kmers[side][kmer]['gap']))
                            if not left in self.half_mers:
                                self.half_mers[left] = {}
                            self.half_mers[left][right] = kmer 
                            left = reverse_complement(left)
                            right = reverse_complement(right)
                            if not right in self.half_mers:
                                self.half_mers[right] = {}
                            self.half_mers[right][left] = kmer
        print(len(self.kmers), 'kmers')
        self.export_accelerator_input()
        self.round_robin()

    def reduce(self):
        c = config.Configuration()
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                line = json_file.readline()
                while line:
                    kmer = line[:line.find(':')]
                    count = int(line[line.find(':') + 1:])
                    self.kmers[kmer]['count'] += count
                    line = json_file.readline()
        print('exporting C kmers...')
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        # output kmers per track
        self.tracks = {}
        for kmer in self.kmers:
            side = self.kmers[kmer]['side']
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner': {}, 'outer': {}}
                self.tracks[track][side][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        # merge
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'gapped_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class GappedKmersIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'GappedKmersIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = CountUniqueGappedKmersJob
    _kmer_type = 'gapped'

    @staticmethod
    def launch(**kwargs):
        job = GappedKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        c = config.Configuration()
        print(cyan(track_name))
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            n = len(kmers['inner']) + len(kmers['outer'])
            for side in ['outer']:
                for kmer in kmers[side]:
                    if not kmer in self.lp_kmers:
                        self.lp_kmers[kmer] = {
                            'gap': kmers[side][kmer]['gap'],
                            'side': side,
                            'type': self._kmer_type,
                            'count': kmers[side][kmer]['count'],
                            'tracks': {},
                            'reference': kmers[side][kmer]['tracks'][track_name],
                            'actual_gap': kmers[side][kmer]['actual_gap'] if 'actual_gap' in kmers[side][kmer] else -1,
                        }
                    self.lp_kmers[kmer]['tracks'][track_name] = 1#kmers[side][kmer]['track']
        path = os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner': {kmer: self.lp_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.lp_kmers, kmers['inner']))},
                    'outer': {kmer: self.lp_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.lp_kmers, kmers['outer']))},
                }, json_file, indent = 4, sort_keys = True)
        return 'gapped_kmers_' + track_name + '.json'

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            kmer['residue'] = 0
            kmer['coverage'] = 42

    def generate_linear_program(self):
        print('generating linear program')
        c = config.Configuration()
        globals()['cplex'] = __import__('cplex')
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        # the coverage of each event
        names = [''] * len(self.tracks)
        for track in self.tracks:
            tokens = track.split('_')
            names[self.tracks[track]['index']] = 'c' + tokens[1]
        problem.variables.add(names = names,
            ub = [1.0] * len(self.tracks),
        )
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            lb = [(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in self.lp_kmers],
            #obj = [1.0] * len(self.lp_kmers),
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
                    rhs = [kmer['count'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['G']
                )
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [-kmer['count'] + sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['G']
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
                    rhs = [kmer['count']],
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
