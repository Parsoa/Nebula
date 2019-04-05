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

from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
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

    def transform(self, track, track_name):
        print(cyan(track_name))
        c = config.Configuration()
        gapped_kmers = self.extract_boundary_gapped_kmers(track)
        path = os.path.join(self.get_current_job_directory(), track_name  + '.json')
        with open(path, 'w') as json_file:
            json.dump(gapped_kmers, json_file, sort_keys = True, indent = 4)
        return track_name + '.json'

    # These kmers ARE NOT CANONICAL
    def extract_boundary_gapped_kmers(self, track):
        c = config.Configuration()
        sequence = track.extract_base_sequence()
        begin = track.slack
        end = len(sequence) - track.slack
        gapped_kmers = {}
        #
        b = track.sequence[begin - c.hsize - 2 - c.ksize: begin + 3 + c.hsize + c.ksize]
        kmer = track.sequence[begin - c.hsize - 2: begin + 3 + c.hsize]
        prefix = track.sequence[begin - c.hsize - 2 - c.ksize: begin - c.hsize - 2]
        suffix = track.sequence[begin + 3 + c.hsize: begin + 3 + c.hsize + c.ksize]
        gapped_kmers[kmer] = {'left': prefix, 'right': suffix, 'side': 'inner'}
        #
        e = track.sequence[end - c.hsize - 2 - c.ksize: end + 3 + c.hsize + c.ksize]
        kmer = track.sequence[end - c.hsize - 2: end + 3 + c.hsize]
        prefix = track.sequence[end - c.hsize - 2 - c.size: end - c.hsize - 2]
        suffix = track.sequence[end + 3 + c.hsize: end + 3 + c.hsize + c.ksize]
        gapped_kmers[kmer] = {'left': prefix, 'right': suffix, 'side': 'inner'}
        #
        kmer = track.sequence[begin - 2 - c.hsize: begin + 3] + track.sequence[end - 2: end + 3 + c.hsize]
        prefix = track.sequence[begin - c.hsize - 2 - c.ksize: begin - c.hsize - 2]
        suffix = track.sequence[end + 3 + c.hsize: end + 3 + c.hsize + c.ksize]
        gapped_kmers[kmer] = {'left': prefix, 'right': suffix, 'side': 'outer'}
        return gapped_kmers

# ============================================================================================================================ #
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
        self.gapped_kmers = {} 
        tracks = self.load_previous_job_results()
        for track in tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                gapped_kmers = json.load(json_file)
                for kmer in gapped_kmers:
                    k = canonicalize(kmer)
                    k = k[:c.hsize] + k[-c.hsize:]
                    if not k in self.gapped_kmers:
                        self.gapped_kmers[k] = gapped_kmers[kmer]
                        self.gapped_kmers[k]['tracks'] = {}
                    if not track in self.gapped_kmers[k]:
                        self.gapped_kmers[k]['tracks'][track] = 0
                    self.gapped_kmers[k]['tracks'][track] += 1
        self.round_robin(tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            _gapped_kmers = json.load(json_file)
        gapped_kmers = {} 
        for kmer in _gapped_kmers:
            if 'N' in kmer:
                continue
            k = canonicalize(kmer)
            k = k[:c.hsize] + k[-c.hsize:]
            gapped_kmers[k] = self.gapped_kmers[k]
        path = os.path.join(self.get_current_job_directory(), track_name  + '.json')
        with open(path, 'w') as json_file:
            json.dump(gapped_kmers, json_file, sort_keys = True, indent = 4)
        return track_name + '.json'

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
                for kmer in kmers:
                    self.kmers[kmer] = kmers[kmer]
                    self.kmers[kmer]['count'] = 0
                    left = kmer[:c.hsize]
                    right = kmer[-c.hsize:]
                    self.half_mers[left] = True
                    self.half_mers[reverse_complement(right)] = True
        print('scoring', len(self.kmers), 'kmers')
        self.chroms = extract_whole_genome()
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
            for k in range(32, 38):
                kmer = canonicalize(sequence[index: index + k])
                kmer = kmer[:c.hsize] + kmer[-c.hsize:]
                if kmer in self.kmers:
                    if self.kmers[kmer]['side'] == 'inner' and k != c.ksize + 5:
                        continue
                    self.kmers[kmer]['count'] += 1

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def merge_counts(self):
        c = config.Configuration()
        print('merging kmer counts ...')
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    self.kmers[kmer]['count'] += batch[kmer]['count']
                    if 'loci' in batch[kmer]:
                        self.kmers[kmer]['loci'].update(batch[kmer]['loci'])

    def reduce(self):
        self.merge_counts()
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)

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
        m = 0
        for track in tracks:
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers:
                    if kmers[kmer]['count'] == 0:
                        left = kmer[:c.hsize]
                        right = kmer[-c.hsize:]
                        self.kmers[kmer] = {'tracks': kmers[kmer]['tracks'], 'count': {}, 'side': kmers[kmer]['side']}
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
        print(len(tracks), 'total events')
        print(len(self.kmers), 'kmers')
        print(m, 'events with kmers')
        self.export_accelerator_input()
        self.round_robin()

    def export_accelerator_input(self):
        with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
            json.dump(self.half_mers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'pre_gapped_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)

    def merge_count(self, kmer, tokens):
        for i in range(0, 11):
            count = tokens[i] 
            self.kmers[kmer]['count'][i] += count

    def reduce(self):
        c = config.Configuration()
        self.tracks = {}
        self.merge_counts()
        print(len(self.kmers))
        keys = list(self.kmers.keys())
        for kmer in keys:
            gap = max(self.kmers[kmer]['count'].items(), key = operator.itemgetter(1))[0]
            total = sum(map(lambda x: self.kmers[kmer]['count'][x], self.kmers[kmer]['count']))
            self.kmers[kmer]['gap'] = gap if self.kmers[kmer]['count'][gap] != 0 else -1
            self.kmers[kmer]['count'] = self.kmers[kmer]['count'][gap]
            self.kmers[kmer]['total'] = total
            self.kmers[kmer]['actual_gap'] = -1
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)
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
                for kmer in kmers:
                    if kmers[kmer]['gap'] != -1:
                        left = kmer[:c.hsize]
                        right = kmer[-c.hsize:]
                        self.kmers[kmer] = kmers[kmer]
                        self.kmers[kmer]['count'] = 0
                        self.kmers[kmer]['doubt'] = 0
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

    def merge_count(self, kmer, tokens):
        count = tokens[0] 
        self.kmers[kmer]['count'] += count

    def reduce(self):
        c = config.Configuration()
        self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        # output kmers per track
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {} 
                self.tracks[track][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)

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
            for kmer in kmers:
                if not kmer in self.lp_kmers:
                    self.lp_kmers[kmer] = {
                        'gap': kmers[kmer]['gap'],
                        'side': kmers[kmer]['side'],
                        'type': self._kmer_type,
                        'count': kmers[kmer]['count'],
                        'tracks': {},
                        'weight': 1.0,
                        'reference': kmers[kmer]['tracks'][track_name],
                        'actual_gap': kmers[kmer]['actual_gap'] if 'actual_gap' in kmers[kmer] else -1,
                    }
                self.lp_kmers[kmer]['tracks'][track_name] = 1
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
            kmer['coverage'] = c.coverage
            kmer['residue'] = 0
            kmer['count'] = min(kmer['count'], kmer['coverage'] * kmer['reference'])

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

    def round_genotype(self, c):
        if c > self.b2:
            return (1.0, '11')
        elif c > self.b1:
            return (0.5, '10')
        else:
            return (0.0, '00')

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
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.reduce)
