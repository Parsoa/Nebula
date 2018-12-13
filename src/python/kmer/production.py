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
    gapped,
    counter,
    reduction,
    simulator,
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

class MixIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'MixIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = MixIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def execute(self):
        c = config.Configuration()
        self.prepare()
        self.create_output_directories()
        self.gapped_kmers = []
        self.gapped_tracks = {}
        self.inner_kmers = []
        self.inner_tracks = {}
        #
        self.load_gapped_kmers()
        self.load_inner_kmers()
        self.merge_kmers()
        self.solve()
        exit()

    def load_gapped_kmers(self):
        print(cyan('=============================================================================================='))
        print(cyan('Loading gapped kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'r') as json_file:
                payload = json.load(json_file)
                self.gapped_kmers = payload['gapped_kmers']
                self.gapped_tracks = payload['tracks']
            return
        self.gapped_kmers_solver = gapped.GappedKmersIntegerProgrammingJob()
        self.gapped_kmers_solver.create_output_directories()
        self.gapped_kmers_solver.load_inputs()
        self.gapped_kmers_solver.distribute_workload()
        self.gapped_kmers_solver.wait_for_children()
        self.gapped_kmers = copy.deepcopy(self.gapped_kmers_solver.index_kmers())
        self.gapped_tracks = copy.deepcopy(self.gapped_kmers_solver.index_tracks())
        with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
            json.dump({'gapped_kmers': self.gapped_kmers, 'tracks': self.gapped_tracks}, json_file, indent = 4)
        self.gapped_kmers_solver.calculate_residual_coverage()
        self.gapped_kmers_solver.solve()
        self.gapped_tracks = copy.deepcopy(self.gapped_kmers_solver.solve())

    def load_inner_kmers(self):
        print(cyan('=============================================================================================='))
        print(cyan('Loading inner kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'r') as json_file:
                payload = json.load(json_file)
                self.inner_kmers = payload['inner_kmers']
                self.inner_tracks = payload['tracks']
            return
        self.inner_kmers_solver = reduction.LociIndicatorKmersIntegerProgrammingJob()
        self.inner_kmers_solver.create_output_directories()
        self.inner_kmers_solver.load_inputs()
        self.inner_kmers_solver.distribute_workload()
        self.inner_kmers_solver.wait_for_children()
        self.inner_kmers = copy.deepcopy(self.inner_kmers_solver.index_kmers())
        self.inner_tracks = copy.deepcopy(self.inner_kmers_solver.index_tracks())
        with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
            json.dump({'inner_kmers': self.inner_kmers, 'tracks': self.inner_tracks}, json_file, indent = 4)
        self.inner_kmers_solver.calculate_residual_coverage()
        self.inner_kmers_solver.solve()

    def merge_kmers(self):
        c = config.Configuration()
        self.tracks = {}
        self.tracks.update(self.inner_tracks)
        self.tracks.update(self.gapped_tracks)
        n = 0
        tmp = sorted([t for t in self.tracks])
        self.lp_kmers = self.gapped_kmers + self.inner_kmers
        for track in tmp:
            self.tracks[track] = {'index': n,
                'inner_kmers': self.inner_tracks[track]['kmers'] if track in self.inner_tracks else [],
                'gapped_kmers': self.gapped_tracks[track]['kmers'] if track in self.gapped_tracks else [],
            }
            n += 1
            with open(os.path.join(self.get_current_job_directory(), 'kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)
        self.lp_kmers = self.inner_kmers + self.gapped_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump({'kmers': self.lp_kmers}, json_file, indent = 4)
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('Solving composite model on'), blue(len(self.lp_kmers)), cyan('kmers and '), blue(len(self.tracks)), cyan('tracks'))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['residue'] = 0 if kmer['type'] == 'gapped' else kmer['reference'] - r
            kmer['coverage'] = c.coverage if kmer['type'] == 'gapped' else kmer['coverage']
            l = len(self.tracks[kmer['tracks'].keys()[0]]['inner_kmers'])
            #kmer['weight'] = 1.0 if kmer['type'] != 'gapped' else 2 * l if l != 0 else 1.0
            l = l if l else 1
            kmer['weight'] = 2 * l if kmer['type'] == 'gapped' else 0.5 + abs(kmer['count'] - (kmer['coverage'] / 2.0)) / (kmer['coverage'] / 2.0) if kmer['reference'] == 1 else 0.5

    def generate_linear_program(self):
        c = config.Configuration()
        import cplex
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
        # the real-valued error parameter for kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            #lb = [-1 * kmer['weight'] * (kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in self.lp_kmers]
            lb = [-100000000 for kmer in self.lp_kmers]
        )
        # absolute value of the kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [1.0 for index, kmer in enumerate(self.lp_kmers)]
        )
        for index, kmer in enumerate(self.lp_kmers):
            self.add_error_absolute_value_constraints(problem, index)
            if kmer['type'] == 'inner':
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
                ind.append(len(self.tracks) + index) # Objective
                val = list(map(lambda track: kmer['weight'] * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
                val.append(1.0) #Objective
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['weight'] * (kmer['count'] - kmer['coverage'] * kmer['residue'])],
                    senses = ['E']
                )
            if kmer['type'] == 'gapped' and kmer['side'] == 'outer':
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: -1 * kmer['weight'] * kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['weight'] * (kmer['count'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))))],
                    senses = ['G']
                )
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: kmer['weight'] * kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['weight'] * (-kmer['count'] + sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))))],
                    senses = ['G']
                )
        return problem

    def add_snp_linear_constraint(self, problem):
        offset = len(self.tracks) + 2 * len(self.lp_kmers)
        n = 0
        for index, kmer in enumerate(self.inner_kmers):
            ref = kmer['reference']
            problem.variables.add(names = ['s' + str(index) + 'L' + str(i) for i in range(0, ref)],
                obj = [1.0] * ref,
                lb  = [0.0] * ref,
                ub  = [1.0] * ref,
            )
            offset += ref
        for track in self.tracks:
            ind = [len(self.tracks) + 2 * len(self.lp_kmers) + index for index in self.tracks[track][kmers]]
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = [1.0] * len(ind),
                )],
                rhs = [math.floor(len(ind) / 3)],
                senses = ['L']
            )
    
    def solve(self):
        programming.IntegerProgrammingJob.solve(self)
        for track in self.tracks:
            self.tracks[track]['kmers'] = {'inner_kmers': {}, 'gapped_kmers': {}}
            self.tracks[track]['errors'] = {'inner_kmers': [], 'gapped_kmers': []}
            self.tracks[track]['indices'] = {'inner_kmers': [], 'gapped_kmers': []}
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['kmers'][kmer['type'] + '_kmers'][kmer['kmer']] = kmer
                self.tracks[track]['errors'][kmer['type'] + '_kmers'].append(abs(self.errors[index]))
                self.tracks[track]['indices'][kmer['type'] + '_kmers'].append(index)
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)
        self.plot_confidence_parameters()

    def plot_confidence_parameters(self):
        x = []
        lengths = []
        errors_std = []
        errors_mean = []
        counts_std = []
        counts_mean = []
        numbers = []
        error_ratios = []
        scores = []
        n = 0
        with open(os.path.join(self.get_current_job_directory(), 'confidence.bed'), 'w') as bed_file:
            for track in self.tracks:
                #if self.tracks[track]['actual_genotype'] == '10' and self.tracks[track]['lp_genotype'] == '11':
                #    continue
                if len(self.tracks[track]['kmers']['inner_kmers']) >= 1:
                    m = 0.01
                    for index in self.tracks[track]['indices']['inner_kmers']:
                        kmer = self.lp_kmers[index]
                        m += max(abs(kmer['count'] - kmer['coverage'] * kmer['residue']),\
                            abs(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])))
                    error_ratios.append(sum(self.tracks[track]['errors']['inner_kmers']) / m)
                    #x.append(self.tracks[track]['actual_genotype'] + '_as_' + self.tracks[track]['lp_genotype'])
                    if self.tracks[track]['actual_genotype'] == '10' and self.tracks[track]['lp_genotype'] == '11':
                        n += 1
                        x.append('True')
                    else:
                        x.append('True' if self.tracks[track]['actual_genotype'] == self.tracks[track]['lp_genotype'] else 'False')
                    t = bed.track_from_name(track)
                    lengths.append(t.end - t.begin)
                    errors_std.append(statistics.std(self.tracks[track]['errors']['inner_kmers']))
                    errors_mean.append(statistics.mean(self.tracks[track]['errors']['inner_kmers']))
                    counts_std.append(statistics.std([self.lp_kmers[index]['count'] for index in self.tracks[track]['indices']['inner_kmers']]))
                    counts_mean.append(statistics.mean([self.lp_kmers[index]['count'] for index in self.tracks[track]['indices']['inner_kmers']]))
                    numbers.append(len(self.tracks[track]['indices']['inner_kmers']))
                    scores.append(self.calculate_confidence_score(self.tracks[track]))
                    bed_file.write(t.chrom + '\t' +
                        str(t.begin) + '\t' +
                        str(t.end) + '\t' +
                        str(x[-1]) + '\t' +
                        str(scores[-1]) + '\n')
        print(n)
        visualizer.violin(x, numbers, 'numbers_distribution', self.get_current_job_directory(), 'Prediction', 'Error')
        visualizer.violin(x, counts_std, 'counts_std_distribution', self.get_current_job_directory(), 'Prediction', 'Error')
        visualizer.violin(x, counts_mean, 'counts_mean_distribution', self.get_current_job_directory(), 'Prediction', 'Error')
        visualizer.violin(x, errors_std, 'errors_std_distribution', self.get_current_job_directory(), 'Prediction', 'Error')
        visualizer.violin(x, errors_mean, 'errors_mean_distribution', self.get_current_job_directory(), 'Prediction', 'Error')
        visualizer.violin(x, error_ratios, 'error_ratio_distribution', self.get_current_job_directory(), 'Prediction', 'Error')
        visualizer.violin(x, lengths, 'event_lenght_distribution', self.get_current_job_directory(), 'Prediction', 'Error')
        visualizer.violin(x, scores, 'confidence_score_distribution', self.get_current_job_directory(), 'Prediction', 'Error')

    def calculate_confidence_score(self, track):
        score = 0
        for error in track['errors']['inner_kmers']:
            score += 1.0 / (1 + error)
        for error in track['errors']['gapped_kmers']:
            score += 10 / (1 + error)
            #score += float(len(track['kmers']['inner_kmers'])) / (1 + error)
        return score

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
