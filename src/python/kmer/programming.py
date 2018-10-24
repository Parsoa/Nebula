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
from kmer.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractInnerKmersJob(map_reduce.Job):

    _name = 'ExtractInnerKmersJob'
    _category = 'programming'
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
        extract_whole_genome()
        self.load_reference_counts_provider()
        self.tracks = self.load_tracks() 
        self.round_robin(self.tracks, filter_func = lambda track: track.end - track.begin > 1000000)

    def transform(self, track, track_name):
        print(cyan(track_name))
        c = config.Configuration()
        inner_kmers = track.extract_inner_kmers(counter = self.reference_counts_provider.get_kmer_count, count = 10, n = 1000, overlap = False, canonical = True)
        kmers = {
            'unique_inner_kmers': {kmer: {'track': inner_kmers[kmer], 'reference': self.reference_counts_provider.get_kmer_count(kmer)} for kmer in list(filter(lambda x: self.reference_counts_provider.get_kmer_count(x) == 1, inner_kmers))},
            'non_unique_inner_kmers': {kmer: {'track': inner_kmers[kmer], 'reference': self.reference_counts_provider.get_kmer_count(kmer)} for kmer in list(filter(lambda x: self.reference_counts_provider.get_kmer_count(x) > 1, inner_kmers))},
        }
        if len(inner_kmers) == 0:
            print(red('skipping', track_name, 'no inner kmers found'))
            return None
        name = 'inner_kmers_' + track_name  + '.json'
        path = os.path.join(self.get_current_job_directory(), name) 
        with open(path, 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4)
        return name

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountInnerKmersJob(counter.BaseExactCountingJob):
    
    _name = 'CountInnerKmersJob'
    _category = 'programming'
    _previous_job = ExtractInnerKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountInnerKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = self.load_previous_job_results()
        for track in self.tracks:
            print(track)
            with open(os.path.join(self.get_previous_job_directory(), self.tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['unique_inner_kmers']:
                    if not kmer in self.kmers:
                        self.kmers[kmer] = {'count': 0, 'track': track, 'reference': kmers['unique_inner_kmers'][kmer]['reference']}
                        self.kmers[reverse_complement(kmer)] = {'count': 0, 'track': track, 'reference': kmers['unique_inner_kmers'][kmer]['reference']}
        self.round_robin()

    def reduce(self):
        self.kmers = self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        self.tracks = {}
        for kmer in self.kmers:
            track = self.kmers[kmer]['track']
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

    _name = 'IntegerProgrammingJob'
    _category = 'programming'
    _previous_job = CountInnerKmersJob
    _kmer_type = 'unique_inner'

    @staticmethod
    def launch(**kwargs):
        job = IntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        tracks = self.load_previous_job_results()
        self.round_robin(tracks)
        self.counts_provider = counttable.DictionaryCountsProvider(json.load(open(os.path.join(self.get_previous_job_directory(), 'kmers.json'))))
        self.lp_kmers = {}

    def transform(self, track, track_name):
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            if len(kmers) == 0:
                print('no unique inner kmers found for', red(track_name))
                return None
            for kmer in kmers:
                count = self.counts_provider.get_kmer_count(str(kmer))
                self.lp_kmers[kmer] = {
                    'type': self._kmer_type,
                    'count': count,
                    'tracks': {track_name: 1},
                    'reference': kmers[kmer]['reference']
                }
            novel_kmers = {}
        path = os.path.join(self.get_current_job_directory(), 'unique_inner_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {kmer: self.lp_kmers[kmer] for kmer in kmers}, json_file, indent = 4, sort_keys = True)
        return path

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.lp_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        c = config.Configuration()
        self.index_kmers()
        self.index_tracks()
        self.calculate_residual_coverage()
        print('exporting kmers...')
        with open(os.path.join(self.get_current_job_directory(), 'lp_kmers.json'), 'w') as json_file:
            json.dump(self.lp_kmers, json_file, indent = 4, sort_keys = True)
        print('generating linear program...')
        self.solve()

    def index_kmers(self):
        c = config.Configuration()
        self.tracks = {}
        self.lp_kmers = []
        index = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if not os.path.isfile(path):
                continue
            #print(path)
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                #print(len(kmers))
                n = 0
                for kmer in kmers:
                    if not kmer in index:
                        index[kmer] = len(self.lp_kmers)
                        self.lp_kmers.append(copy.deepcopy(kmers[kmer]))
                        self.lp_kmers[len(self.lp_kmers) - 1]['kmer'] = kmer
                        for track in kmers[kmer]['tracks']:
                            if not track in self.tracks:
                                self.tracks[track] = {}
                        n += 1
        print(green(len(self.lp_kmers)), 'kmers')
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.lp_kmers, json_file, indent = 4)
        return self.lp_kmers

    def index_tracks(self):
        n = 0
        tmp = sorted([t for t in self.tracks])
        for track in tmp:
            self.tracks[track].update({'index': n, 'kmers': []})
            n += 1
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['kmers'].append(index)
        print(len(self.tracks), 'tracks')
        self.find_tracks_with_no_signal()
        with open(os.path.join(self.get_current_job_directory(), 'tracks.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, indent = 4)
        return self.tracks

    # the portion of a kmer's coverage in reference genome that is outside deletions
    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['residue'] = kmer['reference'] - r
            kmer['coverage'] = c.coverage

    def generate_linear_program(self):
        c = config.Configuration()
        globals()['cplex'] = __import__('cplex')
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        # the coverage of each event
        for track in self.tracks:
            tokens = track.split('_')
            problem.variables.add(names = ['c' + str(tokens[1])],
                ub = [1.0],
                lb = [0.0],
            )
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            ub = [kmer['count'] - kmer['coverage'] * kmer['residue'] for kmer in self.lp_kmers],
            lb = [kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks']) for kmer in self.lp_kmers],
        )
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [1.0] * len(self.lp_kmers),
        )
        #self.add_snp_linear_constraint(problem)
        n = 0
        start = time.time()
        offset = len(self.tracks) + 2 * len(self.lp_kmers)
        for index, kmer in enumerate(self.lp_kmers):
            # TxR + E = C - 
            #ref = kmer['reference']
            ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
            #ind += [offset + i for i in range(0, ref)] #SNP
            ind.append(len(self.tracks) + index) # Objective
            val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
            #val += [-kmer['coverage']] * ref #SNP
            val.append(1.0) #Objective
            #offset += ref
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']],
                senses = ['E']
            )
            #val = [-v for v in val[:-1]]
            #val.append(1.0)
            #problem.linear_constraints.add(
            #    lin_expr = [cplex.SparsePair(
            #        ind = ind,
            #        val = val,
            #    )],
            #    rhs = [-kmer['count'] + kmer['coverage'] * kmer['residue']],
            #    senses = ['E']
            #)
            self.add_error_absolute_value_constraints(problem, index)
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(self.lp_kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
                print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
        return problem

    # We allow up to 3% of the kmers to be affected by SNPs
    def add_snp_linear_constraint(self, problem):
        offset = len(self.tracks) + 2 * len(self.lp_kmers)
        n = 0
        for index, kmer in enumerate(self.lp_kmers):
            ref = kmer['reference']
            problem.variables.add(names = ['s' + str(index) + 'L' + str(i) for i in range(0, ref)],
                obj = [1.0] * ref,
                lb = [0.0] * ref,
                ub = [1.0] * ref,
                #types = [problem.variables.type.integer] * ref,
            )
            offset += ref
        # SNP constraints
        ind = [i for i in range(len(self.tracks) + 2 * len(self.lp_kmers), offset)]
        problem.linear_constraints.add(
            lin_expr = [cplex.SparsePair(
                ind = ind,
                val = [1.0] * len(ind),
            )],
            rhs = [math.ceil(0.03 * len(ind))],
            senses = ['L']
        )

    def add_snp_indicator_variables(self, problem):
        # SNP indicator variables
        problem.variables.add(names = ['i' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            types = [problem.variables.type.binary] * len(self.lp_kmers),
        )
        for track in self.tracks:
            tokens = track.split('_')
            start = int(tokens[1])
            end = int(tokens[2])
            offset = len(self.tracks) + 2 * len(self.lp_kmers)
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = [offset + kmer for kmer in self.tracks[track]['kmers']],
                    val = [1.0] * len(self.tracks[track]['kmers'])
                )],
                rhs = [int(math.floor(0.03 * len(self.tracks[track]['kmers'])))],
                senses = ['L']
            )

    def add_error_absolute_value_constraints(self, problem, index):
        if not 'cplex' in globals():
            globals()['cplex'] = __import__('cplex')
        problem.linear_constraints.add(
            lin_expr = [cplex.SparsePair(
                ind = [len(self.tracks) + len(self.lp_kmers) + index, len(self.tracks) + index],
                val = [1.0, 1.0],
            )],
            rhs = [0],
            senses = ['G']
        )
        problem.linear_constraints.add(
            lin_expr = [cplex.SparsePair(
                ind = [len(self.tracks) + len(self.lp_kmers) + index, len(self.tracks) + index],
                val = [1.0, -1.0],
            )],
            rhs = [0],
            senses = ['G']
        )

    def solve(self):
        problem = self.generate_linear_program()
        problem.write(os.path.join(self.get_current_job_directory(), 'program.lp'))
        problem.solve()
        solution = problem.solution.get_values()
        with open(os.path.join(self.get_current_job_directory(), 'solution.json'), 'w') as json_file:
            json.dump({'variables': problem.solution.get_values()}, json_file, indent = 4, sort_keys = True)
        l = []
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            for track in self.tracks:
                tokens = track.split('_')
                start = int(tokens[1])
                end = int(tokens[2])
                l.append(end - start)
                index = self.tracks[track]['index']
                s = int(round(2 * solution[index]))
                g = '(0, 0)' if s == 2 else '(1, 0)' if s == 1 else '(1, 1)'
                bed_file.write(tokens[0] + '\t' + #0
                            tokens[1] + '\t' + #1
                            tokens[2] + '\t' + #2
                            str(g) + '\t' + #3
                            #str(d) + '\t' + #4
                            str(solution[index]) + '\n') #4
        visualizer.histogram(l, "event lenght", self.get_current_job_directory(), "event lenght", "number of events")

    def find_tracks_with_no_signal(self):
        c = config.Configuration()
        if not c.simulation:
            return
        all_tracks = {}
        n = 0
        for track in self.load_tracks():
            all_tracks[track] = n
            n += 1
        x = []
        with open(os.path.join(self.get_current_job_directory(), 'no_signal.bed'), 'w') as bed_file:
            for track in all_tracks:
                if not track in self.tracks:
                    tokens = track.split('_')
                    x.append(int(tokens[2]) - int(tokens[1]))
                    bed_file.write(tokens[0] + '\t' + #0
                                tokens[1] + '\t' + #1
                                tokens[2] + '\n') #2
        visualizer.histogram(x = x, name = 'events_without_signal_length_distribution', x_label = 'event length', y_label = 'number of events', step = 1, path = self.get_current_job_directory())

    def plot(self, _):
        pass
#
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
