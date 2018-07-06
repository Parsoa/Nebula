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
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import cplex
import numpy
import pybedtools

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
        job = ExtractInnerKmersJob(job_name = 'ExtractInnerKmersJob_', previous_job_name = 'ExtractBreakPointsJob_', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        print('here', self.get_current_job_directory())
        c = config.Configuration()
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.bedtools = {str(track): track for track in pybedtools.BedTool(c.bed_file)}
        print(len(self.bedtools))
        self.round_robin(self.bedtools, lambda track: re.sub(r'\s+', '_', str(track).strip()).strip(), lambda track: track.end - track.start > 1000000) 
        #tracks = self.load_previous_job_results()
        #self.round_robin(tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        tokens = track_name.split('_')
        b = bed.BedTrack(tokens[0], int(tokens[1]), int(tokens[2]))
        sv = self.get_sv_type()(b)
        inner_kmers = sv.get_inner_kmers(self.reference_counts_provider.get_kmer_count, count = 10, n = 1000)
        kmers = {
            'unique_inner_kmers': {kmer: {'track': inner_kmers[kmer], 'count': self.reference_counts_provider.get_kmer_count(kmer)} for kmer in list(filter(lambda x: self.reference_counts_provider.get_kmer_count(x) == 1, inner_kmers))},
            'inner_kmers': {kmer: {'track': inner_kmers[kmer], 'count': self.reference_counts_provider.get_kmer_count(kmer)} for kmer in list(filter(lambda x: self.reference_counts_provider.get_kmer_count(x) > 1, inner_kmers))},
            'novel_kmers': {},
        }
        #with open(track, 'r') as json_file:
        #    break_points = json.load(json_file)
        #    for break_point in break_points:
        #        if break_point.startswith('('):
        #            kmers['novel_kmers'][break_point] = break_points[break_point]['novel_kmers']
        if len(inner_kmers) == 0:
            print(red('skipping', track_name, 'no inner kmers found'))
            return None
        path = os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track_name  + '.json') 
        with open(path, 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class IntegerProgrammingJob(map_reduce.BaseGenotypingJob):

    @staticmethod
    def launch(**kwargs):
        #job = IntegerProgrammingJob(job_name = 'IntegerProgramming_', previous_job_name = 'ExtractInnerKmersJob_', category = 'programming', batch_file_prefix = 'inner_kmers', **kwargs)
        #job.execute()
        job = IntegerProgrammingJob(job_name = 'IntegerProgramming_', previous_job_name = 'ExtractInnerKmersJob_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        tracks = self.load_previous_job_results()
        print(cyan(len(tracks)), 'tracks')
        #items = tracks.items()
        #tracks = {item[0]: item[1] for item in items[0 : 20]}
        self.round_robin(tracks)
        self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.inner_kmers = {}
        self.novel_kmers = {}

    def transform(self, track, track_name):
        with open(track, 'r') as json_file:
            kmers = json.load(json_file)
            #inner_kmers = kmers[self.batch_file_prefix]
            inner_kmers = kmers['unique_inner_kmers']
            if len(inner_kmers) != 0:
                print('no inner kmers found for', red(track_name))
                return None
            inner_kmers = kmers['inner_kmers']
            for kmer in inner_kmers:
                if not kmer in self.inner_kmers:
                    self.inner_kmers[kmer] = {
                        'count': self.counts_provider.get_kmer_count(kmer),
                        'tracks': {},
                        'residue': self.reference_counts_provider.get_kmer_count(kmer)
                    }
                self.inner_kmers[kmer]['tracks'][track_name] = inner_kmers[kmer]['track']
            novel_kmers = {}
            break_points = kmers['novel_kmers']
            for break_point in break_points:
                for kmer in break_points[break_point]:
                    if not kmer in self.novel_kmers:
                        self.novel_kmers[kmer] = {
                            'count': self.counts_provider.get_kmer_count(kmer),
                            'track': track_name,
                            'break_point': break_point,
                        }
                        novel_kmers[kmer] = self.novel_kmers[kmer]
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner_kmers': {kmer: self.inner_kmers[kmer] for kmer in inner_kmers},
                    'novel_kmers': {kmer: self.novel_kmers[kmer] for kmer in novel_kmers},
                }, json_file, indent = 4, sort_keys = True)
        return path

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + str(self.index) + '.json'), 'w')
        json.dump({'inner_kmers': self.inner_kmers, 'novel_kmers': self.novel_kmers}, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        self.inner_kmers = []
        self.novel_kmers = []
        self.tracks = {}
        index = {'inner_kmers': {}, 'novel_kmers': {}}
        inner_kmer_index = 0
        novel_kmer_index = 0
        track_index = 0
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + str(i) + '.json')
            if not os.path.isfile(path):
                continue
            print(path)
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['inner_kmers']:
                    if not kmer in index['inner_kmers']:
                        index['inner_kmers'][kmer] = inner_kmer_index
                        self.inner_kmers.append(kmers['inner_kmers'][kmer])
                        self.inner_kmers[inner_kmer_index]['kmer'] = kmer 
                        inner_kmer_index += 1
                    for track in kmers['inner_kmers'][kmer]['tracks']:
                        self.inner_kmers[index['inner_kmers'][kmer]]['tracks'][track] = kmers['inner_kmers'][kmer]['tracks'][track]
                        if not track in self.tracks:
                            self.tracks[track] = {
                                'index': track_index,
                                'break_points': [],
                                'novel_kmers': []
                            }
                            track_index += 1
                for kmer in kmers['novel_kmers']:
                    if not kmer in index['novel_kmers']:
                        index['novel_kmers'][kmer] = novel_kmer_index
                        self.novel_kmers.append(kmers['novel_kmers'][kmer])
                        self.novel_kmers[novel_kmer_index]['kmer'] = kmer 
                        novel_kmer_index += 1
                        track = kmers['novel_kmers'][kmer]['track']
                        if not track in self.tracks:
                            self.tracks[track] = {
                                'index': track_index,
                                'break_points': [],
                                'novel_kmers': []
                            }
                            track_index += 1
                        self.tracks[track]['break_points'].append({ kmers['novel_kmers'][kmer]['break_point']: 0 })
                        self.tracks[track]['novel_kmers'].append(index['novel_kmers'][kmer])
        print('calculating residual coverage for', green(len(index['inner_kmers'])), 'kmers...')
        for kmer in self.inner_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['residue'] -= r
        print('exporting kmers...')
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_kmers.json'), 'w') as json_file:
            json.dump({'inner_kmers': self.inner_kmers, 'novel_kmers': self.novel_kmers}, json_file, indent = 4, sort_keys = True)
        print('generating linear program...')
        self.solve()

    def generate_linear_program(self):
        c = config.Configuration()
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        self.incorporate_inner_kmers(problem)
        #self.incorporate_novel_kmers(problem, problem.variables.get_num())
        return problem

    #def incorporate_novel_kmers(self, problem, offset):
    #    # the real-valued error parameter for novel_kmers
    #    problem.variables.add(names = ['e' + str(offset + index) for index, kmer in enumerate(self.novel_kmers)],
    #        types = ['C'] * len(self.novel_kmers),
    #        ub = [kmer['count'] for kmer in self.novel_kmers])
    #    b = 0
    #    # break-point selectors
    #    for track in self.tracks:
    #        self.tracks[track]['break_points'] = [{break_point: offset + len(self.novel_kmers) + b + index} for index, break_point in enumerate(self.tracks[track]['break_points'])]
    #        problem.variables.add(names = ['b' + str(offset + len(self.novel_kmers) + b + index) for index, break_point in enumerate(self.tracks[track]['break_points'])],
    #            types = ['B'] * len(self.tracks[track]['break_points']))
    #        b += len(self.tracks[track]['break_points'])
    #    # novel kmer minimization term
    #    for track in self.tracks:
    #        problem.variables.add(names = ['l' + str(offset + len(self.novel_kmers) + b + index) for index, break_point in enumerate(self.tracks[track]['break_points'])],
    #            types = ['C'] * len(self.tracks[track]['break_points']),
    #            obj = [1] * len(self.tracks[track]['break_points']))
    #        b += len(self.tracks[track]['break_points'])
    #        # break point indicators should add up to one
    #        problem.linear_constraints.add(
    #                lin_expr = [cplex.SparsePair(
    #                    ind = [self.tracks[track]['break_points'][break_point] for break_point in self.tracks[track]['break_points']],
    #                    val = [1.0] * len(self.tracks[track]['break_points']
    #                )],
    #                rhs = [1],
    #                senses = ['E']
    #            )
    #        for break_point in self.tracks[track]['break_points']:
    #            problem.quadratic_constraints.add(
    #                quad_expr = [cplex.SparseTriple(
    #                    ind1 = [self.tracks[track]['break_points'][break_point]] * len(self.tracks[track]['novel_kmers']),
    #                    ind2 = [offset + index for index in list(filter(lambda k: self.novel_kmers[k]['break_point'] == break_point, self.tracks[track]['novel_kmers']],
    #                    val  = [-1] * len(self.tracks[track]['novel_kmers']),
    #                ])
    #                lin_expr = [cplex.SparseTuple(
    #                    ind = [problem.variables.get_num() - 1],
    #                    val = [1]
    #                )],
    #                rhs = [0],
    #                senses = ['E']
    #            )
    #    for kmer in self.novel_kmers:
    #        b = self.tracks[kmer['track']]['break_points'][self.tracks[kmer['track']]['break_points'].index(kmer['break_point'])]
    #        problem.quadratic_constraints.add(
    #                lin_expr = [cplex.SparseTuple(
    #                    ind = [b],
    #                    val = [c.coverage]
    #                )]
    #                quad_expr = [cplex.SparseTriple(
    #                    ind1 = [kmer['track']],
    #                    ind2 = [len(self.tracks) + 2 * len(self.inner_kmers) + b],
    #                    val  = [-1 * c.coverage]
    #                )]
    #                rhs = [kmer['count']],
    #                senses = ['E']
    #            )

    def incorporate_inner_kmers(self, problem):
        # the coverage of each event
        problem.variables.add(names = ['c' + str(self.tracks[track]['index']) for track in self.tracks],
            ub = [1.0] * len(self.tracks))
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.inner_kmers)],
            types = ['C'] * len(self.inner_kmers),
            ub = [kmer['count'] for kmer in self.inner_kmers],
            lb = [-1 * c.coverage * sum(kmer['tracks'][track] for track in kmer['tracks']) + kmer['count'] - c.coverage * kmer['residue'] for kmer in self.inner_kmers])
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.inner_kmers)],
            obj = [1.0] * len(self.inner_kmers),
            types = ['C'] * len(self.inner_kmers))
        # constraints
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.inner_kmers):
            ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # C
            ind.append(len(self.tracks) + index) #E + str(index)
            val = list(map(lambda track: c.coverage * kmer['tracks'][track], kmer['tracks']))
            val.append(1.0)
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [kmer['count'] - c.coverage * kmer['residue']],
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
        obj = 0
        for i in range(len(self.tracks), len(self.tracks) + len(self.inner_kmers)):
            obj += solution[i]
        aggregate_count = sum(list(map(lambda kmer: kmer['count'], self.inner_kmers)))
        print('error ratio:', float(obj) / aggregate_count)
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_merge.bed'), 'w') as bed_file:
            for track in self.tracks:
                tokens = track.split('_')
                s = round(2 * solution[self.tracks[track]['index']])
                s = '(0, 0)' if s == 2 else '(1, 0)' if s == 1 else '(1, 1)'
                bed_file.write(tokens[0] + '\t' + tokens[1] + '\t' + tokens[2] + '\t' +  s + '\t' + str(solution[self.tracks[track]['index']]) + '\t' + self.batch_file_prefix + '\n')


    def plot(self, _):
        counts = [kmer['count'] for kmer in self.inner_kmers]
        visualizer.histogram(counts, self.batch_file_prefix, self.get_current_job_directory(), x_label = 'number of times kmer appears in sample', y_label = 'number of kmers')

    def get_previous_job_directory(self):
        return os.path.abspath(os.path.join(map_reduce.Job.get_output_directory(self), self.previous_job_name[:-1]))

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class IntegerProgrammingStatsJob(map_reduce.BaseGenotypingJob):

    @staticmethod
    def launch(**kwargs):
        #job = IntegerProgrammingStatsJob(job_name = 'IntegerProgramming_', previous_job_name = 'IntegerProgramming_', category = 'programming', batch_file_prefix = 'inner_kmers', **kwargs)
        #job.execute()
        job = IntegerProgrammingStatsJob(job_name = 'IntegerProgramming_', previous_job_name = 'IntegerProgramming_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        print(self.batch_file_prefix)
        print(self.get_current_job_directory())
        zygosities = ['00', '10', '11']
        x = []
        y = []
        for z in zygosities:
            for w in zygosities:
                with open(os.path.join(self.get_current_job_directory(), z + '_as_' + w + '.bed')) as bed_file:
                    lines = bed_file.readlines()
                    for line in lines:
                        tokens = line.split()
                        if tokens[6] == self.batch_file_prefix:
                            x.append(z)
                            print(z, w)
                            y.append(float(tokens[5]) + random.randint(1, 10) * 0.001)
        visualizer.violin(x, y, 'LP_' + self.batch_file_prefix, self.get_current_job_directory(), x_label = 'real genotype', y_label = 'LP value')
        exit()

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
