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

from kmer import (
    bed,
    config,
    counter,
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import acora
import cplex
import numpy
import pybedtools

from Bio import pairwise2

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractInnerKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractInnerKmersJob(job_name = 'ExtractInnerKmersJob_', previous_job_name = 'MostLikelyBreakPointsJob_', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.bedtools = {str(track): track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'all.bed'))}
        self.round_robin(self.bedtools, lambda track: re.sub(r'\s+', '_', str(track).strip()).strip(), lambda track: track.end - track.start > 1000000) 

    def transform(self, track, track_name):
        c = config.Configuration()
        sv = self.get_sv_type()(bed.track_from_name(track_name))
        inner_kmers = sv.get_inner_kmers(counter = self.reference_counts_provider.get_kmer_count, count = 10, n = 1000, overlap = False, canonical = True)
        novel_kmers = sv.get_boundary_kmers(begin = 0, end = 0, counter = self.reference_counts_provider.get_kmer_count, count = 1)
        l = len(inner_kmers)
        inner_kmers = {kmer: inner_kmers[kmer] for kmer in filter(lambda k: k not in novel_kmers and reverse_complement(k) not in novel_kmers, inner_kmers)}
        if l != len(inner_kmers):
            print(yellow(track_name))
        kmers = {
            'unique_inner_kmers': {kmer: {'track': inner_kmers[kmer], 'reference': self.reference_counts_provider.get_kmer_count(kmer)} for kmer in list(filter(lambda x: self.reference_counts_provider.get_kmer_count(x) == 1, inner_kmers))},
            'inner_kmers': {kmer: {'track': inner_kmers[kmer], 'reference': self.reference_counts_provider.get_kmer_count(kmer)} for kmer in list(filter(lambda x: self.reference_counts_provider.get_kmer_count(x) > 1, inner_kmers))},
            'novel_kmers': novel_kmers,
        }
        if len(inner_kmers) == 0:
            print(red('skipping', track_name, 'no inner kmers found'))
            return None
        path = os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track_name  + '.json') 
        with open(path, 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountInnerKmersJob(counter.SimulationExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountInnerKmersJob(job_name = 'LocationAwareCountInnerKmersJob_', previous_job_name = 'ExtractInnerKmersJob_', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = self.load_previous_job_results()
        index = {'inner_kmers': {}, 'novel_kmers': {}}
        for track in self.tracks:
            print(track)
            with open(self.tracks[track], 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['inner_kmers']:
                    if not kmer in self.kmers:
                        self.kmers[kmer] = {'count': 0, 'track': track, 'reference': kmers['inner_kmers'][kmer]['reference']}
                        self.kmers[reverse_complement(kmer)] = {'count': 0, 'track': track, 'reference': kmers['inner_kmers'][kmer]['reference']}
                for kmer in kmers['unique_inner_kmers']:
                    if not kmer in self.kmers:
                        self.kmers[kmer] = {'count': 0, 'track': track, 'reference': kmers['inner_kmers'][kmer]['reference']}
                        self.kmers[reverse_complement(kmer)] = {'count': 0, 'track': track, 'reference': kmers['inner_kmers'][kmer]['reference']}
        self.round_robin()

    # delete this once done
    def transform(self, track, track_name):
        c = config.Configuration()
        self.fastq_file = open(track, 'r')
        for read, name in self.parse_fastq():
            kmers = extract_canonical_kmers(read)
            name = name[1:]
            tokens = name.split('_')
            for kmer in kmers:
                if kmer in self.kmers:
                    t = bed.track_from_name(self.kmers[kmer]['track'])
                    if tokens[0] == t.chrom and int(tokens[1]) >= t.start and int(tokens[1]) < t.end:
                        print(self.kmers[kmer]['track'], name)
                        self.kmers[kmer]['count'] += 1

    def reduce(self):
        self.kmers = self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        self.tracks = {}
        for kmer in self.kmers:
            track = self.kmers[kmer]['track']:
            if not track in self.tracks:
                self.tracks[track] = {}
            self.tracks[track][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track + '.json'), 'w') as track_file:
                json.dump(self.tracks[track], track_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'inner_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class IntegerProgrammingJob(map_reduce.BaseGenotypingJob):

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = IntegerProgrammingJob(job_name = 'IntegerProgrammingJob_', previous_job_name = 'CountInnerKmersJob_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', **kwargs)
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
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.inner_kmers = {}
        self.novel_kmers = {}

    def transform(self, track, track_name):
        with open(track, 'r') as json_file:
            kmers = json.load(json_file)
            inner_kmers = {}
            #inner_kmers.update(kmers['inner_kmers'])
            inner_kmers.update(kmers['unique_inner_kmers'])
            if len(inner_kmers) == 0:
                print('no inner kmers found for', red(track_name))
                return None
            for kmer in inner_kmers:
                if not kmer in self.inner_kmers:
                    count = self.counts_provider.get_kmer_count(str(kmer))
                    self.inner_kmers[kmer] = {
                        'type': 'unique_inner',
                        'count': count,
                        'tracks': {},
                        'reference': self.reference_counts_provider.get_kmer_count(kmer)
                    }
                if kmer in self.inner_kmers:
                    self.inner_kmers[kmer]['tracks'][track_name] = inner_kmers[kmer]['track']
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'unique_inner_kmers': {kmer: self.inner_kmers[kmer] for kmer in inner_kmers},
                }, json_file, indent = 4, sort_keys = True)
        return path

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + str(self.index) + '.json'), 'w')
        json.dump({'unique_inner_kmers': self.inner_kmers}, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        c = config.Configuration()
        self.index_kmers()
        self.index_tracks()
        self.calculate_residual_coverage()
        print('exporting kmers...')
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_kmers.json'), 'w') as json_file:
            json.dump({'unique_inner_kmers': self.inner_kmers}, json_file, indent = 4, sort_keys = True)
        print('generating linear program...')
        self.solve()

    def index_kmers(self):
        c = config.Configuration()
        self.tracks = {}
        self.inner_kmers = []
        index = {'inner_kmers': {}}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + str(i) + '.json')
            if not os.path.isfile(path):
                continue
            print(path)
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['inner_kmers']:
                    if not kmer in index['inner_kmers']:
                        index['inner_kmers'][kmer] = len(self.inner_kmers)
                        self.inner_kmers.append(copy.deepcopy(kmers['inner_kmers'][kmer]))
                        self.inner_kmers[len(self.inner_kmers) - 1]['kmer'] = kmer
                        for track in kmers['inner_kmers'][kmer]['tracks']:
                            if not track in self.tracks:
                                self.tracks[track]['unique_inner_kmers'] = 0
                            self.tracks[track]['unique_inner_kmers'] += kmers['inner_kmers'][kmer]['tracks']
        print(green(len(self.inner_kmers)), 'inner kmers')
        return self.inner_kmers

    def index_tracks(self):
        n = 0
        tmp = sorted([t for t in self.tracks])
        for track in tmp:
            self.tracks[track].update({'index': n})
            n += 1
        print(len(self.tracks), 'tracks')
        self.find_tracks_with_no_signal()
        return self.tracks

    # the portion of a kmer's coverage in reference genome that is outside deletions
    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.inner_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['residue'] = kmer['reference'] - r
            kmer['coverage'] = c.coverage

    def generate_linear_program(self):
        c = config.Configuration()
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        self.incorporate_inner_kmers(problem)
        #self.incorporate_novel_kmers(problem, problem.variables.get_num())
        return problem

    def incorporate_inner_kmers(self, problem):
        # the coverage of each event
        for track in self.tracks:
            tokens = track.split('_')
            problem.variables.add(names = ['c' + str(tokens[1])],
                ub = [1.0],
            )
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.inner_kmers)],
            ub = [kmer['count'] - kmer['coverage'] * kmer['residue'] for kmer in self.inner_kmers],
            lb = [kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks']) for kmer in self.inner_kmers]
        )
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.inner_kmers)],
            obj = [1.0] * len(self.inner_kmers),
        )
        #problem.objective.set_quadratic([0.0] * len(self.tracks) + [1.0] * len(self.inner_kmers) + [0.0] * len(self.inner_kmers))
        #problem.objective.set_quadratic([cplex.SparsePair(ind = [index], val = [1.0]) for index, kmer in enumerate(self.inner_kmers)])
        # constraints
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.inner_kmers):
            ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # C
            ind.append(len(self.tracks) + index) #E + str(index)
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
                    ind = [len(self.tracks) + len(self.inner_kmers) + index, len(self.tracks) + index],
                    val = [1.0, 1.0],
                )],
                rhs = [0],
                senses = ['G']
            )
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = [len(self.tracks) + len(self.inner_kmers) + index, len(self.tracks) + index],
                    val = [1.0, -1.0],
                )],
                rhs = [0],
                senses = ['G']
            )
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(self.inner_kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
                print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
        return problem

    def solve(self):
        problem = self.generate_linear_program()
        problem.write(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_program.lp'))
        problem.solve()
        solution = problem.solution.get_values()
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_solution.json'), 'w') as json_file:
            json.dump({'variables': problem.solution.get_values()}, json_file, indent = 4, sort_keys = True)
        #obj = 0
        #for i in range(len(self.tracks), len(self.tracks) + len(self.inner_kmers)):
        #    obj += abs(solution[i])
        #max_error = sum(list(map(lambda kmer: max(abs(kmer['count'] - kmer['coverage'] * kmer['residue']), abs(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks']))), self.inner_kmers)))
        #print('error ratio:', float(obj) / max_error)
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_merge.bed'), 'w') as bed_file:
            for track in self.tracks:
                tokens = track.split('_')
                index = self.tracks[track]['index']
                s = int(round(2 * solution[index]))
                s = '(0, 0)' if s == 2 else '(1, 0)' if s == 1 else '(1, 1)'
                bed_file.write(tokens[0] + '\t' + #0
                            tokens[1] + '\t' + #1
                            tokens[2] + '\t' + #2
                            s + '\t' + #3
                            str(solution[index]) + '\t' + #4
                            #str(len(track['inner_kmers'])) + '\t' + #5
                            self.batch_file_prefix + '\n') #6

    def find_tracks_with_no_signal(self):
        bedtools = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'all.bed'))], key = lambda track: track.start)
        all_tracks = {}
        n = 0
        for b in bedtools:
            track = re.sub(r'\s+', '_', str(b).strip()).strip()
            print(track)
            all_tracks[track] = n
            n += 1
        with open(os.path.join(self.get_current_job_directory(), 'no_signal.bed'), 'w') as bed_file:
            for track in all_tracks:
                if not track in self.tracks:
                    tokens = track.split('_')
                    bed_file.write(tokens[0] + '\t' + #0
                                tokens[1] + '\t' + #1
                                tokens[2] + '\n') #2

    def plot(self, _):
        counts = [kmer['count'] for kmer in self.inner_kmers]
        visualizer.histogram(counts, self.batch_file_prefix, self.get_current_job_directory(), x_label = 'number of times kmer appears in sample', y_label = 'number of kmers')

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        d = 'simulation' if c.simulation else self.category
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../' + d + '/' + bed_file_name + '/' + str(c.ksize) + '/', self.previous_job_name[:-1]))

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #
#
#class IntegerProgrammingStatsJob(map_reduce.BaseGenotypingJob):
#
#    @staticmethod
#    def launch(**kwargs):
#        #job = IntegerProgrammingStatsJob(job_name = 'IntegerProgramming_', previous_job_name = 'IntegerProgramming_', category = 'programming', batch_file_prefix = 'inner_kmers', **kwargs)
#        #job.execute()
#        job = IntegerProgrammingStatsJob(job_name = 'IntegerProgramming_', previous_job_name = 'IntegerProgramming_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', **kwargs)
#        job.execute()
#
#    # ============================================================================================================================ #
#    # MapReduce overrides
#    # ============================================================================================================================ #
#
#    def load_inputs(self):
#        print(self.batch_file_prefix)
#        print(self.get_current_job_directory())
#        self.zygosities = ['00', '10', '11']
#        self.tracks = {}
#        x = []
#        y = []
#        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as merge_file:
#            for z in self.zygosities:
#                self.tracks[z] = {}
#                for w in self.zygosities:
#                    self.tracks[z][w] = []
#                    with open(os.path.join(self.get_current_job_directory(), z + '_as_' + w + '.bed')) as bed_file:
#                        lines = bed_file.readlines()
#                        for line in lines:
#                            line = line.strip()
#                            tokens = line.strip().split('\t')
#                            name = tokens[0] + '_' + tokens[1] + '_' + tokens[2]
#                            print(line)
#                            with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + name + '.json'), 'r') as json_file:
#                                kmers = json.load(json_file)['inner_kmers']
#                                line += '\t'
#                                line += '(0, 0) 'if z == '00' else '(1, 0)' if z == '10' else '(1, 1)'
#                                for kmer in kmers:
#                                    line += '\t' + str(kmers[kmer]['count'])
#                                line += '\n'
#                                merge_file.write(line)
#                            self.tracks[z][w].append(tokens)
#                            print(tokens)
#                            if tokens[6] == self.batch_file_prefix:
#                                x.append(z)
#                                y.append(float(tokens[4]) + random.randint(1, 10) * 0.001)
#        visualizer.violin(x, y, 'LP_' + self.batch_file_prefix, self.get_current_job_directory(), x_label = 'real genotype', y_label = 'LP value')
#        #self.calc_sv_kmer_count_variance_and_genotype_correlation()
#        exit()
#
#    def profile_lp_value_per_number_of_kmers(self):
#        x = []
#        y = []
#        for z in self.zygosities:
#            for track in self.tracks['10'][z]:
#                t = int(track[5])
#                x.append(min((t / 10) * 10, 50))
#                x.append(min((t / 10) * 10, 50))
#                y.append(float(track[4]) + random.randint(1, 10) * 0.001)
#                y.append(float(track[4]) + random.randint(1, 10) * 0.001)
#        print(x)
#        print(y)
#        visualizer.violin(x, y, '10_prediction_per_number_of_kmers', self.get_current_job_directory(), x_label = 'number of kmers', y_label = 'LP value')
#
#    def calc_sv_kmer_count_variance_and_genotype_correlation(self):
#        from scipy.stats.stats import pearsonr
#        x = []
#        y = []
#        for z in self.zygosities:
#            for w in self.zygosities:
#                for track in self.tracks[z][w]:
#                    with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track[0] + '_' + track[1] + '_' + track[2] + '.json')) as json_file:
#                        kmers = json.load(json_file)['inner_kmers']
#                        counts = [kmers[kmer]['count'] for kmer in kmers]
#                        std = statistics.variance(counts)
#                        if z == w:
#                            x.append(1)
#                        else:
#                            x.append(0)
#                        y.append(min(std, 2000))
#        for i, j in zip(x, y):
#            print(i, j)
#        visualizer.violin(x, y, 'variance_per_genotype_correctness', self.get_current_job_directory(), x_label = 'genotype correctness (1 = correct, 0 = wrong)', y_label = 'variance')
#        print(pearsonr(x, y))
#
## ============================================================================================================================ #
## ============================================================================================================================ #
## ============================================================================================================================ #
## ============================================================================================================================ #
## ============================================================================================================================ #
#
#class IterativeIntegerProgrammingJob(IntegerProgrammingJob):
#
#    @staticmethod
#    def launch(**kwargs):
#        c = config.Configuration()
#        job = IterativeIntegerProgrammingJob(job_name = 'IterativeIntegerProgramming_', previous_job_name = 'CountInnerKmersJob_' if c.simulation else 'ExtractInnerKmersJob_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', **kwargs)
#        job.execute()
#
#    def load_inputs(self):
#        c = config.Configuration()
#        zygosities = ['00', '10', '11']
#        self.excluded_tracks = {}
#        for z in zygosities:
#            for w in zygosities:
#                with open(os.path.join(self.get_current_job_directory(), z + '_as_' + w + '.bed')) as bed_file:
#                    lines = bed_file.readlines()
#                    for line in lines:
#                        tokens = line.split()
#                        name = tokens[0] + '_' + tokens[1] + '_' + tokens[2]
#                        if z == '11':
#                            self.excluded_tracks[name] = 0.0
#                        if z == '10':
#                            self.excluded_tracks[name] = 0.5
#                        if z == '00':
#                            self.excluded_tracks[name] = 1.0
#                        #r = random.randint(1, 4)
#                        #if r == 3:
#                        #    if z == '11':
#                        #        self.excluded_tracks[name] = 0.0#1.0 if z == '00' else 0.5 if z == '10' else 0.0
#                        #    if z == '10':
#                        #        self.excluded_tracks[name] = 0.5#1.0 if z == '00' else 0.5 if z == '10' else 0.0
#                        #    if z == '00':
#                        #        self.excluded_tracks[name] = 1.0#1.0 if z == '00' else 0.5 if z == '10' else 0.0
#        tracks = self.load_previous_job_results()
#        #for track in self.excluded_tracks:
#        #    tracks.pop(track, None)
#        self.round_robin(tracks)
#        #self.dic_counts_provider = counttable.DictionaryCountsProvider(json.load(open(os.path.join(self.get_previous_job_directory(), 'kmers.json'))))
#        if c.simulation:
#            self.counts_provider = counttable.DictionaryCountsProvider(json.load(open(os.path.join(self.get_previous_job_directory(), 'kmers.json'))))
#        else:
#            self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
#        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
#        self.inner_kmers = {}
#        self.novel_kmers = {}
#
#    def incorporate_inner_kmers(self, problem):
#        IntegerProgrammingJob.incorporate_inner_kmers(self, problem)
#        for index, track in enumerate(self.tracks):
#            if track['track'] in self.excluded_tracks:
#                problem.variables.set_lower_bounds(index, self.excluded_tracks[track['track']] - 0.01)
#                problem.variables.set_upper_bounds(index, self.excluded_tracks[track['track']] + 0.01)
#
## ============================================================================================================================ #
## ============================================================================================================================ #
## ============================================================================================================================ #
## ============================================================================================================================ #
## ============================================================================================================================ #
#
#class PerEventIntegerProgrammingJob(IntegerProgrammingJob):
#
#    @staticmethod
#    def launch(**kwargs):
#        c = config.Configuration()
#        job = PerEventIntegerProgrammingJob(job_name = 'PerEventIntegerProgrammingJob_', previous_job_name = 'CountInnerKmersJob_' if c.simulation else 'ExtractInnerKmersJob_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', **kwargs)
#        job.execute()
#
#    def incorporate_inner_kmers(self, problem, track):
#        # the coverage of each event
#        tokens = track['track'].split('_')
#        problem.variables.add(names = ['c' + str(tokens[1])], ub = [1.0])
#        # the real-valued error parameter for inner_kmer
#        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(track['inner_kmers'])],
#            types = ['C'] * len(track['inner_kmers']),
#            ub = [kmer['count'] - kmer['coverage'] * kmer['residue'] for index, kmer in enumerate(track['inner_kmers'])],
#            lb = [kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks']) for index, kmer in enumerate(track['inner_kmers'])]
#        )
#        # absolute value of the inner_kmer error parameter
#        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(track['inner_kmers'])],
#            obj = [1.0] * len(track['inner_kmers']),
#            types = ['C'] * len(track['inner_kmers']))
#        # constraints
#        #json_print(track)
#        #print(len(track['inner_kmers']))
#        for index, kmer in enumerate(track['inner_kmers']):
#            #print(kmer)
#            ind = [0]
#            ind.append(1 + index)#E
#            #print(ind)
#            val = [kmer['coverage'] * kmer['tracks'][track['index']]]
#            val.append(1.0)
#            #print(val)
#            problem.linear_constraints.add(
#                lin_expr = [cplex.SparsePair(
#                    ind = ind,
#                    val = val,
#                )],
#                rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']],
#                senses = ['E']
#            )
#            problem.linear_constraints.add(
#                lin_expr = [cplex.SparsePair(
#                    ind = [1 + len(track['inner_kmers']) + index, 1 + index],
#                    val = [1.0, 1.0],
#                )],
#                rhs = [0],
#                senses = ['G']
#            )
#            problem.linear_constraints.add(
#                lin_expr = [cplex.SparsePair(
#                    ind = [1 + len(track['inner_kmers']) + index, 1 + index],
#                    val = [1.0, -1.0],
#                )],
#                rhs = [0],
#                senses = ['G']
#            )
#        return problem
#
#    def solve(self):
#        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_merge.bed'), 'w') as bed_file:
#            for track in self.tracks:
#                print('solving for', track['track'])
#                problem = cplex.Cplex()
#                problem.objective.set_sense(problem.objective.sense.minimize)
#                self.incorporate_inner_kmers(problem, track)
#                problem.write(os.path.join(self.get_current_job_directory(), str(track['track']) + '_program.lp'))
#                problem.solve()
#                #print('exporting lp')
#                solution = problem.solution.get_values()
#                tokens = track['track'].split('_')
#                s = int(round(2 * solution[0]))
#                s = '(0, 0)' if s == 2 else '(1, 0)' if s == 1 else '(1, 1)'
#                bed_file.write(tokens[0] + '\t' + #0
#                            tokens[1] + '\t' + #1
#                            tokens[2] + '\t' + #2
#                            s + '\t' + #3
#                            str(solution[0]) + '\t' + #4
#                            str(len(track['inner_kmers'])) + '\t' + #5
#                            self.batch_file_prefix + '\n') #6
#
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
