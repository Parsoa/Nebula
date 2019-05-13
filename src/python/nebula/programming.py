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

from nebula import (
    bed,
    config,
    counter,
    map_reduce,
    statistics,
    visualizer,
)

from nebula.kmers import *
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

from pulp import *

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class IntegerProgrammingJob(map_reduce.BaseGenotypingJob):

    _name = 'IntegerProgrammingJob'
    _category = 'preprocessing'
    _previous_job = None
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
        self.lp_kmers = {}

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
        with open(os.path.join(self.get_current_job_directory(), 'tracks.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, indent = 4)
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
                debug_log('batch not found:', path)
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers:
                    if not kmer in index:
                        index[kmer] = len(self.lp_kmers)
                        self.lp_kmers.append(copy.deepcopy(kmers[kmer]))
                        self.lp_kmers[len(self.lp_kmers) - 1]['kmer'] = kmer
                        for track in kmers[kmer]['tracks']:
                            if not track in self.tracks:
                                self.tracks[track] = {}
        print(green(len(self.lp_kmers)), 'kmers')
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
        return self.tracks

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            # put an upperbound on a kmer's impact on LP score
            kmer['count'] = min(kmer['count'], kmer['coverage'] * kmer['reference'])
            kmer['residue'] = kmer['reference'] - r

    def add_error_absolute_value_constraints(self, problem, index):
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

    def add_mps_error_absolute_value_constraints(self, problem, variables, index):
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], 1.0), (variables[len(self.tracks) + index], 1.0)])
        problem += LpConstraint(expr, LpConstraintGE, 'abs_1_' + str(index), 0) 
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], 1.0), (variables[len(self.tracks) + index], -1.0)])
        problem += LpConstraint(expr, LpConstraintGE, 'abs_2_' + str(index), 0) 

    def solve(self):
        c = config.Configuration()
        if c.solver == 'coin':
            problem, variables = self.generate_mps_linear_program()
            problem.writeLP(os.path.join(self.get_current_job_directory(), 'program_coin.lp'))
            problem.writeMPS(os.path.join(self.get_current_job_directory(), 'program_coin.mps'))
            command = '/share/hormozdiarilab/Codes/NebulousSerendipity/coin/build/bin/clp ' + os.path.join(self.get_current_job_directory(), 'program_coin.mps') + ' -dualsimplex -solution ' + os.path.join(self.get_current_job_directory(), 'solution.mps')
            output = subprocess.call(command, shell = True)
            self.import_lp_values()
            with open(os.path.join(self.get_current_job_directory(), 'solution_coin.json'), 'w') as json_file:
                json.dump({'variables': self.solution}, json_file, indent = 4, sort_keys = True)
        if c.solver == 'cplex':
            problem = self.generate_linear_program()
            problem.solve()
            problem.write(os.path.join(self.get_current_job_directory(), 'program_cplex.lp'))
            self.solution = problem.solution.get_values()
            with open(os.path.join(self.get_current_job_directory(), 'solution_cplex.json'), 'w') as json_file:
                json.dump({'variables': self.solution}, json_file, indent = 4, sort_keys = True)
        self.export_solution()
        self.verify_genotypes()
        return self.tracks

    # COIN doesn't supply values for certain variables
    def import_lp_values(self):
        c = config.Configuration()
        self.solution = [0.0] * (len(self.tracks) + 2 * len(self.lp_kmers))
        var_index = {}
        regex = re.compile('[^a-zA-Z0-9]')
        for track in self.tracks:
            name = regex.sub('_', track)
            var_index[name] = self.tracks[track]['index']
        with open(os.path.join(self.get_current_job_directory(), 'solution.mps'), 'r') as f:
            status = f.readline()
            objective = f.readline()
            line = f.readline()
            while(line):
                tokens = line.split()
                name = tokens[1]
                index = int(tokens[0])
                value = float(tokens[2])
                if name[0] == 'c':
                    self.solution[var_index[name[1:]]] = value
                else:
                    self.solution[index] = value
                line = f.readline()

    def export_solution(self):
        c = config.Configuration()
        self.errors = self.solution[len(self.tracks):]
        for track in self.tracks:
            index = self.tracks[track]['index']
            self.tracks[track]['coverage'] = self.solution[index]
            self.tracks[track]['lp_kmers'] = []
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['lp_kmers'].append(kmer)
        self.find_rounding_break_points()
        print('Rounding', len(self.tracks), 'tracks')
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            bed_file.write('CHROM\tBEGIN\tEND\tLP_GENOTYPE\tLP_VALUE\tID\n')
            for track in self.tracks:
                t = c.tracks[track]
                index = self.tracks[track]['index']
                g = self.round_genotype(self.solution[index], t.svtype)
                bed_file.write(t.chrom + '\t' +
                    str(t.begin) + '\t' +
                    str(t.end) + '\t' +
                    str(g[1]) + '\t' +
                    str(self.solution[index]) + '\t' + 
                    str(t.id) + '\n')

    def find_rounding_break_points(self):
        c = config.Configuration()
        d_0 = statistics.NormalDistribution(0, min(c.std / 4.0, 3))
        d_1 = statistics.NormalDistribution(c.coverage / 2, c.std / 2)
        d_2 = statistics.NormalDistribution(c.coverage, c.std)
        distributions = [{'inner': d_0, 'gapped': d_2, 'junction': d_2}, {'inner': d_1, 'gapped': d_1, 'junction': d_1}, {'inner': d_2, 'gapped': d_0, 'junction': d_0}]
        m = None
        step = 0.05
        self.b2 = 0.75
        self.b1 = 0.25
        if not c.rounding:
            return
        for b1 in range(int(0.1 / step), int(0.5 / step)):
            for b2 in range(int(0.55 / step), int(1.0 / step)):
                print(b1, b2)
                likelihood = 0
                b_1 = b1 * step
                b_2 = b2 * step
                for track in self.tracks:
                    if self.tracks[track]['coverage'] <= b_1:
                        likelihood += self.estimate_likelihood(track, distributions[0])
                    elif self.tracks[track]['coverage'] > b_1 and self.tracks[track]['coverage'] <= b_2:
                        likelihood += self.estimate_likelihood(track, distributions[1])
                    else:
                        likelihood += self.estimate_likelihood(track, distributions[2])
                if m == None:
                    self.b1 = b_1
                    self.b2 = b_2
                    m = likelihood
                elif likelihood > m:
                    self.b1 = b_1
                    self.b2 = b_2
                    m = likelihood
        print('[0,', cyan(self.b1), ',', yellow(self.b2), ', 1.0]')

    def estimate_likelihood(self, track, distribution):
        c = config.Configuration()
        likelihood = 0.0
        for index, kmer in enumerate(self.tracks[track]['lp_kmers']):
            l = distribution[kmer['type']].log_pmf(kmer['count'] - kmer['residue'] * kmer['coverage'])
            kmer['likelihood'] = l
            likelihood += l
        return likelihood / (index + 1)

    def round_genotype(self, c, svtype):
        if c > self.b2:
            if svtype == 'DEL':
                return (1.0, '00')
            return (1.0, '11')
        elif c > self.b1:
            return (0.5, '10')
        else:
            if svtype == 'DEL':
                return (0.0, '11')
            return (0.0, '00')

    def verify_genotypes(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        output = subprocess.call('verify_genotypes', shell = True, stderr = subprocess.STDOUT, cwd = self.get_current_job_directory())
        for r in ['00', '10', '11']:
            for p in ['00', '10', '11']:
                for track in bed.load_tracks_from_file(os.path.join(self.get_current_job_directory(), r + '_as_' + p + '.bed'), [('lp_genotype', None, str), ('lp_value', None, float)]):
                    g = self.round_genotype(track.lp_value)
                    self.tracks[str(track)]['lp_value'] = track.lp_value
                    self.tracks[str(track)]['lp_rounding'] = g[0]
                    self.tracks[str(track)]['lp_genotype'] = g[1]
                    self.tracks[str(track)]['actual_genotype'] = r
        with open(os.path.join(self.get_current_job_directory(), 'tracks.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, indent = 4)
        #self.plot_lp_values()
        #self.calculate_length_accuracy_correlation()

    def plot_lp_values(self):
        print(len(self.tracks))
        x = []
        for track in self.tracks:
            x.append(self.tracks[track]['lp_value'])
        visualizer.histogram(x, 'lp_value_distribution', self.get_current_job_directory(), 'Value', 'Frequency', 0.05)

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.reduce)
