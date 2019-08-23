from __future__ import print_function

import io
import gc
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

from nebula import (
    bed,
    cgc,
    config,
    counter,
    junction,
    reduction,
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
    preprocessor,
)

from nebula.debug import *
from nebula.kmers import *
from nebula.logger import *
from nebula.chromosomes import *
print = pretty_print

import numpy as np

from pulp import *
from scipy import stats

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UnifiedGenotypingJob(cgc.CgcIntegerProgrammingJob):

    _name = 'UnifiedGenotypingJob'
    _category = 'genotyping'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_genotyping_paths(self):
        #self.paths = ['/share/hormozdiarilab/Codes/NebulousSerendipity/output/clustering/CutlessCandyChr17/' + str(i) for i in range(self.begin, self.end)]
        base = '/share/hormozdiarilab/Codes/NebulousSerendipity/output/cgc/genotyping/july-19-19'
        for (dirpath, dirnames, filenames) in os.walk(base):
            self.paths = [os.path.join(dirpath, dirname) for dirname in dirnames]
            break
        self.paths = sorted(self.paths)[self.begin: self.end]
        print(self.paths)
        self.gold = []
        for path in self.paths:
            base, name = os.path.split(path)
            self.gold.append(os.path.join('/share/hormozdiarilab/Codes/NebulousSerendipity/data/Audano/Simons', 'Audano.SGDP_' + name + '.hg38.ALL.bed'))

    def load_inputs(self):
        self.load_genotyping_paths()
        tracks = json.load(open(os.path.join(self.paths[0], 'CgcCounterJob', 'batch_merge.json' )))
        self.round_robin(tracks)
        self.lp_kmers = {}

    def get_previous_job_directory(self):
        return os.path.join(self.paths[0], 'CgcCounterJob')

    def get_current_job_directory(self):
        return os.path.abspath(os.path.join(self.get_output_directory(), 'UnifiedGenotypingOrchestrator', str(self.genotyping_batch)))

    def transform(self, path, track_name):
        kmers = json.load(open(os.path.join(self.get_previous_job_directory(), path), 'r'))
        c = config.Configuration()
        t = c.tracks[track_name]
        lp_kmers = {}
        for kmer in kmers['inner_kmers']:
            if len(kmers['inner_kmers'][kmer]['loci']) > 3:
                continue
            if any(['chr17' in track for track in kmers['inner_kmers'][kmer]['tracks']]):
                lp_kmers[kmer] = True
                self.lp_kmers[kmer] = kmers['inner_kmers'][kmer]
                self.lp_kmers[kmer]['reduction'] = kmers['inner_kmers'][kmer]['reference']
                self.lp_kmers[kmer]['reference'] = len(kmers['inner_kmers'][kmer]['loci'])
        for kmer in kmers['junction_kmers']:
            if kmers['junction_kmers'][kmer]['source'] == 'assembly':
                continue
            if len(kmers['junction_kmers'][kmer]['loci']) > 1:
                continue
            if any(['chr17' in track for track in kmers['junction_kmers'][kmer]['tracks']]):
                lp_kmers[kmer] = True
                self.lp_kmers[kmer] = kmers['junction_kmers'][kmer]
                self.lp_kmers[kmer]['reduction'] = kmers['junction_kmers'][kmer]['reference']
                self.lp_kmers[kmer]['reference'] = len(kmers['junction_kmers'][kmer]['loci'])
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4)
        return path

    def reduce(self):
        c = config.Configuration()
        self.index_kmers()
        self.index_tracks()
        self.calculate_residual_coverage()
        self.update_counts()
        self.solve()
        self.gather_genotype_statistics()

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
                kmer['svtype'] = c.tracks[track].svtype
            kmer['count'] = []
            kmer['residue'] = kmer['reference'] - r
            kmer['coverage'] = [40] * len(self.paths)
        self.calculate_error_weights()

    def update_counts(self):
        for i, path in enumerate(self.paths):
            print(green(path))
            with open(os.path.join(path, 'CgcCounterJob', 'kmers.json'), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in self.lp_kmers:
                    if any(['chr17' in track for track in kmer['tracks']]):
                        k = kmer['kmer']
                        for kmer_type in ['inner_kmers', 'junction_kmers']:
                            if k in kmers[kmer_type]:
                                kmer['count'].append(min(kmers[kmer_type][k]['count'], kmer['coverage'][i] * kmer['reference']))
            gc.collect()

    def generate_mps_linear_program(self):
        c = config.Configuration()
        problem = LpProblem("Nebula", LpMinimize)
        variables = [None] * len(self.paths) * len(self.tracks)
        n = len(self.paths)
        regex = re.compile('[^a-zA-Z0-9]')
        for p, path in enumerate(self.paths):
            for track in self.tracks:
                name = regex.sub('_', track)
                prefix = 'lp' + str(p) + 'c' + name
                variables[p * len(self.tracks) + self.tracks[track]['index']] = LpVariable(prefix, 0, 1)
        self.add_mps_error_absolute_value_constraints(problem, variables)
        self.add_mps_coverage_diff_absolute_value_constraints(problem, variables)
        expr = LpAffineExpression([(variables[n * (len(self.tracks) + len(self.lp_kmers)) + i], 0 * self.lp_kmers[i % len(self.lp_kmers)]['weight']) for i in range(0, n * len(self.lp_kmers))] + [(variables[j], 1.0 / n) for j in range(n * (len(self.tracks) + 2 * len(self.lp_kmers)), len(variables))])
        problem += expr
        self.add_mps_optimization_constraints(problem, variables)
        return problem, variables

    def add_mps_error_absolute_value_constraints(self, problem, variables):
        # error variables
        for p, path in enumerate(self.paths):
            for index, kmer in enumerate(self.lp_kmers):
                variables.append(LpVariable('p' + str(p) + 'e' + str(index)))
        # absolute value of the error variables
        for p, path in enumerate(self.paths):
            for index, kmer in enumerate(self.lp_kmers):
                variables.append(LpVariable('p' + str(p) + 'l' + str(index)))
        n = 0
        for p, path in enumerate(self.paths):
            for i, kmer in enumerate(self.lp_kmers):
                offset = len(self.paths) * (len(self.tracks) + len(self.lp_kmers))
                expr = LpAffineExpression([(variables[len(self.paths) * (len(self.tracks) + len(self.lp_kmers)) + p * len(self.lp_kmers) + i], 1.0), (variables[len(self.paths) * len(self.tracks) + p * len(self.lp_kmers) + i], 1.0)])
                problem += LpConstraint(expr, LpConstraintGE, 'abs_' + str(p) + '_' + str(n) + '_1', 0)
                expr = LpAffineExpression([(variables[len(self.paths) * (len(self.tracks) + len(self.lp_kmers)) + p * len(self.lp_kmers) + i], 1.0), (variables[len(self.paths) * len(self.tracks) + p * len(self.lp_kmers) + i], -1.0)])
                problem += LpConstraint(expr, LpConstraintGE, 'abs_' + str(p) + '_' + str(n) + '_2', 0)
                n += 1

    def add_mps_coverage_diff_absolute_value_constraints(self, problem, variables):
        offset = len(self.paths) * (len(self.tracks) +  2 * len(self.lp_kmers))
        i = 0
        n = len(self.paths)
        for p in range(0, n):
            for q in range(p + 1, n):
                for k, kmer in enumerate(self.lp_kmers):
                    variables.append(LpVariable('d' + str(i)))
                    expr = LpAffineExpression([(variables[offset + i], 1.0)] + [(variables[p * len(self.tracks) + self.tracks[track]['index']], kmer['coverage'][p]) for track in kmer['tracks']] + [(variables[q * len(self.tracks) + self.tracks[track]['index']], - kmer['coverage'][q]) for track in kmer['tracks']])
                    problem += LpConstraint(expr, LpConstraintGE, '_'.join(['diff', str(k), str(p), str(q), '1']), kmer['count'][q] - kmer['count'][p] if kmer['svtype'] == 'DEL' and kmer['type'] == 'inner' else kmer['count'][p] - kmer['count'][q])
                    expr = LpAffineExpression([(variables[offset + i], 1.0)] + [(variables[q * len(self.tracks) + self.tracks[track]['index']], kmer['coverage'][q]) for track in kmer['tracks']] + [(variables[p * len(self.tracks) + self.tracks[track]['index']], - kmer['coverage'][p]) for track in kmer['tracks']])
                    problem += LpConstraint(expr, LpConstraintGE, '_'.join(['diff', str(k), str(p), str(q), '2']), kmer['count'][p] - kmer['count'][q] if kmer['svtype'] == 'DEL' and kmer['type'] == 'inner' else kmer['count'][q] - kmer['count'][p])
                    i += 1

    def add_mps_optimization_constraints(self, problem, variables):
        c = config.Configuration()
        for p, path in enumerate(self.paths):
            for i, kmer in enumerate(self.lp_kmers):
                indices = list(map(lambda track: p * len(self.tracks) + self.tracks[track]['index'], kmer['tracks']))
                indices.append(len(self.paths) * len(self.tracks) + p * len(self.lp_kmers) + i)
                if kmer['svtype'] == 'INS' or kmer['type'] == 'junction':
                    coeffs = list(map(lambda track: kmer['coverage'][p] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                    coeffs.append(1.0)
                    rhs = kmer['count'][p] - kmer['coverage'][p] * kmer['residue']
                if kmer['svtype'] == 'DEL' and kmer['type'] == 'inner':
                    coeffs = list(map(lambda track: -1 * kmer['coverage'][p] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                    coeffs.append(1.0)
                    rhs = kmer['count'][p] - kmer['coverage'][p] * kmer['residue'] - sum(map(lambda track: kmer['coverage'][p] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(p) + '_' + str(i), rhs)

    def solve(self):
        c = config.Configuration()
        problem, variables = self.generate_mps_linear_program()
        problem.writeLP(os.path.join(self.get_current_job_directory(), 'program_coin.lp'))
        problem.writeMPS(os.path.join(self.get_current_job_directory(), 'program_coin.mps'))
        command = '/share/hormozdiarilab/Codes/NebulousSerendipity/coin/build/bin/clp ' + os.path.join(self.get_current_job_directory(), 'program_coin.mps') + ' -dualsimplex -solution ' + os.path.join(self.get_current_job_directory(), 'solution.mps')
        output = subprocess.call(command, shell = True)
        self.import_lp_values()
        with open(os.path.join(self.get_current_job_directory(), 'solution_coin.json'), 'w') as json_file:
            json.dump({'variables': self.solution}, json_file, indent = 4, sort_keys = True)
        self.export_solution()
        self.export_kmers()
        return self.tracks

    def import_lp_values(self, path = 'solution.mps'):
        c = config.Configuration()
        self.solution = []
        for path in self.paths:
            self.solution.append([0.0] * len(self.tracks))
        var_index = {}
        regex = re.compile('[^a-zA-Z0-9]')
        for p, path in enumerate(self.paths):
            for track in self.tracks:
                name = regex.sub('_', track)
                var_index['lp' + str(p) + 'c' + name] = self.tracks[track]['index']
        with open(os.path.join(self.get_current_job_directory(), 'solution.mps'), 'r') as f:
            status = f.readline()
            objective = f.readline()
            line = f.readline()
            while(line):
                tokens = line.split()
                name = tokens[1]
                value = float(tokens[2])
                if name.startswith('lp'):
                    index = var_index[name]
                    path = int(name[2: name.find('c')])
                    self.solution[path][index] = value
                line = f.readline()

    def export_solution(self):
        c = config.Configuration()
        for track in self.tracks:
            index = self.tracks[track]['index']
            self.tracks[track]['coverage'] = []
            self.tracks[track]['lp_kmers'] = []
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['lp_kmers'].append(kmer)
        for p, path in enumerate(self.paths):
            print('Genotyping sample', path)
            self.find_rounding_break_points()
            with open(os.path.join(self.get_current_job_directory(), 'merge_' + str(p) + '.bed'), 'w') as bed_file:
                with open(os.path.join(path, 'CgcIntegerProgrammingJob', 'union.bed'), 'w') as cluster_file:
                    bed_file.write('#CHROM\tBEGIN\tEND\tLP_GENOTYPE\tLP_VALUE\tID\n')
                    cluster_file.write('#CHROM\tBEGIN\tEND\tLP_GENOTYPE\tLP_VALUE\tID\n')
                    for track in self.tracks:
                        t = c.tracks[track]
                        index = self.tracks[track]['index']
                        g = self.round_genotype(self.solution[p][index], t.svtype)
                        line = t.chrom + '\t' + str(t.begin) + '\t' + str(t.end) + '\t' +\
                            str(g[1]) + '\t' + str(self.solution[p][index]) + '\t' + str(t.id) + '\n'
                        bed_file.write(line)
                        cluster_file.write(line)

    def gather_genotype_statistics(self):
        self.diff = {}
        self.stats = {}
        self.tracks = {}
        self.lp_values = {}
        for name in ['merge', 'union']:
            self.diff[name] = {}
            self.stats[name] = {'genotype': {}, 'likelihood': {}, 'state': {}}
        for p in ['00', '10', '11']:
            for q in ['00', '10', '11']:
                for name in ['merge', 'union']:
                    s = p + '_as_' + q
                    self.stats[name]['genotype'][s] = 0
        for p in ['t', 'f']:
            for q in ['p', 'n']:
                for name in ['merge', 'union']:
                    self.diff[name][p + q] = []
                    self.stats[name]['state'][p + q] = 0
        for i in range(len(self.paths)):
            cwd = os.path.join(self.paths[i], 'CgcIntegerProgrammingJob')
            if i == 0:
                tracks = bed.load_tracks_from_file_as_dict(os.path.join(cwd, 'merge.bed'))
                for track in tracks:
                    self.tracks[track] = {}
                    self.lp_values[track] = []
                    for name in ['merge', 'union']:
                        self.tracks[track][name] = {}
                        for p in ['00', '10', '11']:
                            for q in ['00', '10', '11']:
                                self.tracks[track][name][p + '_as_' + q] = [0, []]
            self.tabulate('merge', cwd, i)
            self.tabulate('union', cwd, i)
        self.diff['diff'] = {'merge': {}, 'union': {}}
        for p in ['t', 'f']:
            for q in ['p', 'n']:
                for name in ['merge', 'union']:
                    self.diff[name][p + q] = set(self.diff[name][p + q])
                self.diff['diff']['merge'][p + q] = list(self.diff['merge'][p + q].difference(self.diff['union'][p + q]))
                self.diff['diff']['union'][p + q] = list(self.diff['union'][p + q].difference(self.diff['merge'][p + q]))
        with open(os.path.join(self.get_current_job_directory(), 'statistics.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'comparison.json'), 'w') as json_file:
            json.dump(self.stats, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'diff.json'), 'w') as json_file:
            json.dump(self.diff['diff'], json_file, indent = 4)

    def tabulate(self, name, cwd, i):
        p = subprocess.Popen(['/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/tabulate.sh', name + '.bed'], cwd = cwd)
        p.wait()
        p = subprocess.Popen(['/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/verify_sim.sh', self.gold[i]], cwd = cwd)
        p.wait()
        for p in ['00', '10', '11']:
            for q in ['00', '10', '11']:
                s = p + '_as_' + q
                tracks = bed.load_tracks_from_file_as_dict(os.path.join(cwd, s + '.bed'))
                for track in tracks:
                    if not 'chr17' in track:
                        continue
                    self.tracks[track][name][s][0] += 1
                    self.tracks[track][name][s][1].append(i)
                    self.stats[name]['genotype'][s] += 1
                    if name == 'merge':
                        self.lp_values[track].append(float(tracks[track].lp_value))
                    if p == q:
                        if '1' in p:
                            self.diff[name]['tp'].append((track, i))
                            self.stats[name]['state']['tp'] += 1
                        else:
                            self.diff[name]['tn'].append((track, i))
                            self.stats[name]['state']['tn'] += 1
                    else:
                        if not '1' in p:
                            self.diff[name]['fp'].append((track, i))
                            self.stats[name]['state']['fp'] += 1
                        else:
                            if '1' in q:
                                self.diff[name]['tp'].append((track, i))
                                self.stats[name]['state']['tp'] += 1
                            else:
                                self.diff[name]['fn'].append((track, i))
                                self.stats[name]['state']['fn'] += 1

    def plot_coverage_count_difference(self, track):
        print(track)
        x = []
        for kmer in self.lp_kmers:
            if track in kmer['tracks']:
                for i in range(self.begin - 1000, self.end - 1000):
                    for j in range(self.begin - 1000, self.end - 1000):
                        if i != j:
                            x.append(abs((self.lp_values[track][i] - self.lp_values[track][j]) * 40 - kmer['count'][i] + kmer['count'][j]))
        if len(x) > 10000:
            x = [x[i] for i in sorted(random.sample(xrange(len(x)), 10000))] 
        visualizer.histogram(x, 'coverage_count_diff_' + track, self.get_current_job_directory(), 'diff', 'count')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UnifiedGenotypingOrchestrator(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'UnifiedGenotypingOrchestrator'
    _category = 'genotyping'
    _previous_job = None 
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = UnifiedGenotypingOrchestrator(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        batch_size = 10
        tracks = {}
        for i in range(0, 400 / batch_size):
            print(yellow('============================================================================'))
            print(yellow('============================================================================'))
            print(yellow('=================================' + str(i) + '======================================='))
            print(yellow('============================================================================'))
            print(yellow('============================================================================'))
            job = UnifiedGenotypingJob(begin = 1000 + i * batch_size, end = 1000 + (i + 1) * batch_size, genotyping_batch = i)
            job.execute()
        exit()
    
    def transform(self, job, index):
        job.execute()
        return None

    def reduce(self):
        pass
