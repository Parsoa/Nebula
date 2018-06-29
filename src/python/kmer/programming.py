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
        job = ExtractInnerKmersJob(job_name = 'ExtractInnerKmersJob_', previous_job_name = '', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        print('here', self.get_current_job_directory())
        c = config.Configuration()
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.bedtools = {str(track): track for track in pybedtools.BedTool(c.bed_file)}
        self.round_robin(self.bedtools, lambda track: re.sub(r'\s+', '_', str(track).strip()).strip(), lambda track: track.end - track.start > 1000000) 

    def transform(self, track, track_name):
        sv = self.get_sv_type()(track)
        c = config.Configuration()
        inner_kmers = sv.get_inner_kmers(self.reference_counts_provider.get_kmer_count, count = 1, n = 1000)
        if len(inner_kmers) == 0:
            return None
        path = os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump({'inner_kmers': inner_kmers}, json_file, sort_keys = True, indent = 4)
        json_file.close()
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
        job = IntegerProgrammingJob(job_name = 'IntegerProgramming_', previous_job_name = 'ExtractInnerKmersJob_', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        tracks = self.load_previous_job_results()
        items = tracks.items()
        #tracks = {item[0]: item[1] for item in items[0 : 20]}
        self.round_robin(tracks)
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.kmers = {}

    def transform(self, track, track_name):
        with open(track, 'r') as json_file:
            inner_kmers = json.load(json_file)['inner_kmers']
            for kmer in inner_kmers:
                if not kmer in self.kmers:
                    self.kmers[kmer] = {
                        'count': self.counts_provider.get_kmer_count(kmer),
                        'tracks': {},
                        'residue': self.reference_counts_provider.get_kmer_count(kmer)
                    }
                self.kmers[kmer]['tracks'][track_name] = inner_kmers[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track_name + '.json'), 'w') as json_file:
            json.dump({kmer: self.kmers[kmer] for kmer in inner_kmers}, json_file, indent = 4, sort_keys = True)
        return True

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        self.kmers = []
        self.tracks = {}
        k_index = 0
        t_index = 0
        k = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + str(i) + '.json')
            print(path)
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers:
                    if not kmer in k:
                        k[kmer] = k_index
                        self.kmers.append(kmers[kmer])
                        self.kmers[k_index]['kmer'] = kmer 
                        k_index += 1
                    for track in kmers[kmer]['tracks']:
                        self.kmers[k[kmer]]['tracks'][track] = kmers[kmer]['tracks'][track]
                        if not track in self.tracks:
                            self.tracks[track] = t_index
                            t_index += 1
        print('calculating residual coverage for', green(len(self.kmers)), '...')
        for kmer in self.kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['residue'] -= r
        print('exporting kmers')
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        print('generating linear program...')
        self.solve()

    def solve(self):
        problem = self.generate_linear_program()
        problem.write(os.path.join(self.get_current_job_directory(), 'program.lp'))
        problem.solve()
        solution = problem.solution.get_values()
        with open(os.path.join(self.get_current_job_directory(), 'solution.json'), 'w') as json_file:
            json.dump({'variables': problem.solution.get_values()}, json_file, indent = 4, sort_keys = True)
        obj = 0
        for i in range(len(self.tracks), len(self.tracks) + len(self.kmers)):
            obj += solution[i]
        aggregate_count = sum(list(map(lambda kmer: kmer['count'], self.kmers)))
        print('error ratio:', float(obj) / aggregate_count)
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            for track in self.tracks:
                tokens = track.split('_')
                s = round(2 * solution[self.tracks[track]])
                s = '(0, 0)' if s == 2 else '(1, 0)' if s == 1 else '(1, 1)'
                bed_file.write(tokens[0] + '\t' + tokens[1] + '\t' + tokens[2] + '\t' +  s + '\n')

    def generate_linear_program(self):
        c = config.Configuration()
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        # the coverage of each event
        problem.variables.add(names = ['c' + str(self.tracks[track]) for track in self.tracks], ub = [1.0] * len(self.tracks))
        # the real-valued error parameter
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.kmers)], types = ['C'] * len(self.kmers),
            ub = [kmer['count'] for kmer in self.kmers],
            lb = [min(-1 * c.coverage * sum([item[1] for item in kmer['tracks'].items()]), kmer['count'] - c.coverage * kmer['residue']) for kmer in self.kmers])

        # the kmer count slack
        #problem.variables.add(names = ['s' + str(index) for index, kmer in enumerate(self.kmers)], ub = [1.0] * len(self.kmers))
        # absolute value of the error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.kmers)], obj = [1.0] * len(self.kmers), types = ['C'] * len(self.kmers))
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.kmers):
            ind = list(map(lambda track: self.tracks[track], kmer['tracks'])) # C
            ind.append(len(self.tracks) + index) #E + str(index)
            #ind.append(len(self.tracks) + len(self.kmers) + index)
            val = list(map(lambda track: c.coverage * kmer['tracks'][track], kmer['tracks']))
            val.append(1.0)
            #val.append(c.coverage * kmer['residue'])
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
                print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
        return problem

    def plot(self, _):
        counts = [kmer['count'] for kmer in self.kmers]
        visualizer.histogram(counts, 'unique_inner_kmers', self.get_current_job_directory(), x_label = 'number of times kmer appears in sample',y_label = 'number of kmers')

    def get_previous_job_directory(self):
        return os.path.abspath(os.path.join(map_reduce.Job.get_output_directory(self), self.previous_job_name[:-1]))

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
