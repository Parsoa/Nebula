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

from nebula import (
    bed,
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
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

import numpy as np

from pulp import *
from scipy import stats

#import plotly.offline as plotly
#import plotly.graph_objs as graph_objs

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcCounterJob(counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'CgcCounterJob'
    _category = 'genotyping'
    _previous_job = preprocessor.MixKmersJob
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = CgcCounterJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.load_kmers()
        self.export_counter_input()
        self.round_robin()

    def load_kmers(self):
        c = config.Configuration()
        path = c.kmers[0]
        with open(path, 'r') as json_file:
            self.kmers = json.load(json_file)
        self.gc_kmers = {}#self.kmers['gc_kmers']
        self.depth_kmers = self.kmers['depth_kmers']
        self.inner_kmers = self.kmers['inner_kmers']
        self.junction_kmers = self.kmers['junction_kmers']
        for path in c.kmers[1:]:
            with open(path, 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['junction_kmers']:
                    if not kmer in self.junction_kmers:
                        self.junction_kmers[kmer] = kmers['junction_kmers'][kmer]
                    else:
                        self.junction_kmers[kmer]['loci'].update(kmers['junction_kmers'][kmer]['loci'])
                        self.junction_kmers[kmer]['tracks'].update(kmers['junction_kmers'][kmer]['tracks'])
                for kmer in kmers['inner_kmers']:
                    if not kmer in self.inner_kmers:
                        self.inner_kmers[kmer] = kmers['inner_kmers'][kmer]
                    else:
                        self.inner_kmers[kmer]['loci'].update(kmers['inner_kmers'][kmer]['loci'])
                        self.inner_kmers[kmer]['tracks'].update(kmers['inner_kmers'][kmer]['tracks'])
        print('Counting...')
        print('Junction kmers:', blue(len(self.junction_kmers)))
        print('Inner kmers:', blue(len(self.inner_kmers)))
        print('Depth kmers:', blue(len(self.depth_kmers)))
        print('GC kmers:', blue(len(self.gc_kmers)))

    def export_counter_input(self):
        if self.resume_from_reduce:
            return
        _kmers = {}
        for kmers in [self.inner_kmers, self.junction_kmers]:
            for kmer in kmers:
                _kmers[kmer] = {}
                _kmers[kmer]['loci'] = {}
                _kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                for locus in kmers[kmer]['loci']:
                    _kmers[kmer]['loci'][locus] = {}
                    _kmers[kmer]['loci'][locus]['masks'] = kmers[kmer]['loci'][locus]['masks']
        for kmer in self.depth_kmers:
            _kmers[kmer] = {}
            _kmers[kmer]['loci'] = {}
        for kmer in self.gc_kmers:
            _kmers[kmer] = {}
            _kmers[kmer]['loci'] = {}
        print('Dumping', green(len(_kmers)), 'kmers for the counter')
        with open(os.path.join(self.get_current_job_directory(), 'pre_inner_kmers.json'), 'w') as json_file:
            json.dump(_kmers, json_file, indent = 4)

    def merge_count(self, kmer, tokens):
        count = tokens[0] 
        total = tokens[1] 
        canon = canonicalize(kmer)
        if canon in self.junction_kmers:
            self.junction_kmers[canon]['count'] += count / 2
            self.junction_kmers[canon]['total'] += total / 2
        if canon in self.inner_kmers:
            self.inner_kmers[canon]['count'] += count / 2
            self.inner_kmers[canon]['total'] += total / 2
        if canon in self.depth_kmers:
            self.depth_kmers[canon]['count'] += count / 2
        if canon in self.gc_kmers:
            self.depth_kmers[canon]['count'] += count / 2

    def estimate_depth_of_coverage(self):
        c = config.Configuration()
        counts = list(map(lambda kmer: self.depth_kmers[kmer]['count'], self.depth_kmers))
        mean = np.mean(counts)
        std = np.std(counts)
        print(len(counts))
        print('mean:', mean)
        print('std:', std)
        # filter outliers
        counts = list(filter(lambda x: x < 3 * mean, counts))
        mean = np.mean(counts)
        std = np.std(counts)
        print(len(counts))
        print('mean:', mean)
        print('std:', std)
        # filter outliers
        counts = list(filter(lambda x: x < 2 * mean, counts))
        mean = np.mean(counts)
        std = np.std(counts)
        print(len(counts))
        print('mean:', mean)
        print('std:', std)
        #
        stats = {'coverage': mean, 'std': std}
        self.estimate_gc_content_coverage()
        with open(os.path.join(self.get_current_job_directory(), 'stats_' + str(c.ksize) + '.json'), 'w') as json_file:
            json.dump(stats, json_file, sort_keys = True, indent = 4)
        return stats

    def estimate_gc_content_coverage(self):
        c = config.Configuration()

    def reduce(self):
        c = config.Configuration()
        self.merge_counts()
        self.export_counted_kmers()
        self.tracks = {}
        print('Exporting genotyping tracks...')
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.junction_kmers:
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        #for track in self.tracks:
        #    with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
        #        json.dump(self.tracks[track], json_file, indent = 4)
        print(len(self.tracks), 'tracks')
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)
        job = ExportGenotypingKmersJob()
        job.tracks = self.tracks
        job.execute()
        return self.estimate_depth_of_coverage()

    def export_counted_kmers(self):
        c = config.Configuration()
        print('Exporting counted kmers...')
        #with open(os.path.join(self.get_current_job_directory(), 'gc_kmers'), 'w') as json_file:
        #    json.dump(self.gc_kmers, json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'depth_kmers'), 'w') as json_file:
        #    json.dump(self.depth_kmers, json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
        #    json.dump(self.inner_kmers, json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'junction_kmers.json'), 'w') as json_file:
        #    json.dump(self.junction_kmers, json_file, indent = 4)
        kmers = {}
        kmers['gc_kmers'] = self.gc_kmers
        kmers['depth_kmers'] = self.depth_kmers
        kmers['inner_kmers'] = self.inner_kmers
        kmers['junction_kmers'] = self.junction_kmers
        if c.cgc:
            name = 'kmers_' + (c.fastq.split('/')[-1] if c.fastq else c.bam.split('/')[-1]) + '.json'
            path = os.path.join(os.getcwd(), name)
        else:
            name = 'kmers.json'
            path = os.path.join(self.get_current_job_directory(), name)
        with open(path, 'w') as json_file:
            json.dump(kmers, json_file, indent = 4)


# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'CgcIntegerProgrammingJob'
    _category = 'genotyping'
    _previous_job = CgcCounterJob

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

    def load_inputs(self):
        c = config.Configuration()
        tracks = self.load_previous_job_results()
        self.round_robin(tracks)
        self.resume_from_reduce = False
        self.lp_kmers = {}

    def transform(self, path, track_name):
        #print(green(track_name))
        kmers = json.load(open(os.path.join(self.get_previous_job_directory(), path), 'r'))
        c = config.Configuration()
        t = c.tracks[track_name]
        lp_kmers = {}
        for kmer in kmers['inner_kmers']:
            if len(kmers['inner_kmers'][kmer]['loci']) > 3:
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['inner_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['inner_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['inner_kmers'][kmer]['loci'])
        for kmer in kmers['junction_kmers']:
            if kmers['junction_kmers'][kmer]['source'] == 'assembly':
                continue
            if len(kmers['junction_kmers'][kmer]['loci']) > 1:
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['junction_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['junction_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['junction_kmers'][kmer]['loci'])
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4)
        return path

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
                kmer['svtype'] = c.tracks[track].svtype
            kmer['coverage'] = c.coverage
            kmer['residue'] = kmer['reference'] - r
            kmer['weight'] = 1.0
            kmer['count'] = min(kmer['count'], kmer['coverage'] * kmer['reference'])
        self.calculate_error_weights()

    def calculate_error_weights(self):
        for track in self.tracks:
            self.add_weights_to_track(track)

    def add_weights_to_track(self, track):
        l = self.tracks[track]['kmers']
        inner = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'inner']
        junction = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'junction']
        for kmers in [inner, junction]:
            for kmer in kmers:
                kmer['weight'] = float(len(inner) + len(junction)) / (2.0 * len(kmers))

    def generate_linear_program(self):
        c = config.Configuration()
        globals()['cplex'] = __import__('cplex')
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        # the coverage of each event
        names = [''] * len(self.tracks)
        regex = re.compile('[^a-zA-Z0-9]')
        for track in self.tracks:
            name = regex.sub('', track)
            names[self.tracks[track]['index']] = 'c' + name
        problem.variables.add(names = names,
            ub = [1.0] * len(self.tracks),
        )
        # the real-valued error parameter for kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            lb = [-100000000 for kmer in self.lp_kmers]
        )
        # absolute value of the kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [kmer['weight'] for index, kmer in enumerate(self.lp_kmers)]
        )
        print('adding constraints')
        for index, kmer in enumerate(self.lp_kmers):
            self.add_error_absolute_value_constraints(problem, index)
            if kmer['type'] == 'junction':
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
                ind.append(len(self.tracks) + index) # Objective
                val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
                val.append(1.0) #Objective
                rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']]
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']],
                    senses = ['E']
                )
            if kmer['type'] == 'inner':
                for track in kmer['tracks']:
                    if c.tracks[track].svtype == 'INS':
                        ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
                        ind.append(len(self.tracks) + index) # Objective
                        val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
                        val.append(1.0) #Objective
                        rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']]
                        problem.linear_constraints.add(
                            lin_expr = [cplex.SparsePair(
                                ind = ind,
                                val = val,
                            )],
                            rhs = rhs,
                            senses = ['E']
                        )
                    else:
                        ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
                        ind.append(len(self.tracks) + index) # Objective
                        val = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
                        val.append(1.0) #Objective
                        rhs = [kmer['count'] - kmer['coverage'] * kmer['residue'] - sum(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))]
                        problem.linear_constraints.add(
                            lin_expr = [cplex.SparsePair(
                                ind = ind,
                                val = val,
                            )],
                            rhs = rhs,
                            senses = ['E']
                        )
                    break
        print('completed program generation')
        return problem

    def generate_mps_linear_program(self):
        c = config.Configuration()
        problem = LpProblem("Nebula", LpMinimize)
        i = 0
        names = [''] * len(self.tracks)
        variables = [None] * (len(self.tracks) + 2 * len(self.lp_kmers))
        regex = re.compile('[^a-zA-Z0-9]')
        for track in self.tracks:
            name = regex.sub('_', track)
            variables[self.tracks[track]['index']] = LpVariable('c' + name, 0, 1)
            problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintLE, 'c' + name + '_ub', 1.0)
            problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintGE, 'c' + name + '_lb', 0.0)
            i += 1
        # error variables
        for index, kmer in enumerate(self.lp_kmers):
            variables[i] = LpVariable('e' + str(index))
            i += 1
        # absolute value of the error variables
        for index, kmer in enumerate(self.lp_kmers):
            variables[i] = LpVariable('l' + str(index))
            i += 1
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], kmer['weight']) for index, kmer in enumerate(self.lp_kmers)])
        problem += expr
        for i, kmer in enumerate(self.lp_kmers):
            self.add_mps_error_absolute_value_constraints(problem, variables, i)
            indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
            indices.append(len(self.tracks) + i)
            if kmer['svtype'] == 'INS' or kmer['type'] == 'junction':
                coeffs = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['count'] - kmer['coverage'] * kmer['residue']
            if kmer['svtype'] == 'DEL' and kmer['type'] == 'inner':
                coeffs = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['count'] - kmer['coverage'] * kmer['residue'] - sum(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
            expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
            problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
        return problem, variables

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcInnerKmersIntegerProgrammingJob(CgcIntegerProgrammingJob):

    _name = 'CgcInnerKmersIntegerProgrammingJob'
    _category = 'genotyping'
    _previous_job = CgcCounterJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcInnerKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    def transform(self, path, track_name):
        print(green(track_name))
        kmers = json.load(open(os.path.join(self.get_previous_job_directory(), path), 'r'))
        c = config.Configuration()
        lp_kmers = {}
        for kmer in kmers['inner_kmers']:
            #if len(kmers['junction_kmers']) != 0:
            #    break
            if len(kmers['inner_kmers'][kmer]['loci']) > 3:
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['inner_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['inner_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['inner_kmers'][kmer]['loci'])
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        print(path)
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4, sort_keys = True)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcJunctionKmersIntegerProgrammingJob(CgcIntegerProgrammingJob):

    _name = 'CgcJunctionKmersIntegerProgrammingJob'
    _category = 'preprocessing'
    _previous_job = CgcCounterJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcJunctionKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    def transform(self, path, track_name):
        print(green(track_name))
        kmers = json.load(open(os.path.join(self.get_previous_job_directory(), path), 'r'))
        c = config.Configuration()
        lp_kmers = {}
        for kmer in kmers['junction_kmers']:
            #if len(kmers['inner_kmers']) != 0:
            #    break
            if kmers['junction_kmers'][kmer]['source'] == 'assembly':
                continue
            if len(kmers['junction_kmers'][kmer]['loci']) > 1:
                continue
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['junction_kmers'][kmer]
            self.lp_kmers[kmer]['reduction'] = kmers['junction_kmers'][kmer]['reference']
            self.lp_kmers[kmer]['reference'] = len(kmers['junction_kmers'][kmer]['loci'])
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        print(path)
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4, sort_keys = True)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcCoverageCorrectingIntegerProgrammingJob(CgcIntegerProgrammingJob):

    _name = 'CgcCoverageCorrectingIntegerProgrammingJob'
    _category = 'genotyping'
    _previous_job = CgcCounterJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcCoverageCorrectingIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def index_tracks(self):
        n = 0
        tmp = sorted([t for t in self.tracks])
        for track in tmp:
            self.tracks[track].update({'index': n, 'kmers': [], 'coverage': [], 'genotypes': []})
            n += 1
        for index, kmer in enumerate(self.lp_kmers):
            kmer['original_count'] = kmer['count']
            for track in kmer['tracks']:
                self.tracks[track]['kmers'].append(index)
        print(len(self.tracks), 'tracks')
        return self.tracks

    def reduce(self):
        c = config.Configuration()
        self.index_kmers()
        self.index_tracks()
        helper = CgcCoverageCorrectingIntegerProgrammingJob.LpHelper()
        helper.tracks = self.tracks
        helper.lp_kmers = self.lp_kmers
        helper.num_iterations = 10
        helper.execute()

    # ============================================================================================================================ #
    # Helper for running several LPs
    # ============================================================================================================================ #

    class LpHelper(CgcIntegerProgrammingJob):

        _name = 'CgcCoverageCorrectingIntegerProgrammingJob'
        _category = 'genotyping'
        _previous_job = CgcCounterJob

        def load_inputs(self):
            iterations = {}
            self.resume_from_reduce = False
            for iteration in range(self.num_iterations):
                iterations[iteration] = True
            self.round_robin(iterations)

        def transform(self, track, track_name):
            print(green(track_name))
            self.calculate_residual_coverage()
            self.solve()
            exit()

        def calculate_residual_coverage(self):
            c = config.Configuration()
            for kmer in self.lp_kmers:
                r = 0
                for track in kmer['tracks']:
                    r += kmer['tracks'][track]
                kmer['coverage'] = {track: c.coverage for track in kmer['tracks']} if self.index == 0 else {track: np.random.normal(c.coverage, c.std) for track in kmer['tracks']}
                kmer['residue'] = c.coverage * (kmer['reference'] - r) if self.index == 0 else sum([np.random.normal(c.coverage, c.std) for i in range(kmer['reference'] - r)])
                kmer['weight'] = 1.0
                kmer['count'] = min(kmer['original_count'], sum([kmer['coverage'][track] * kmer['tracks'][track] for track in kmer['tracks']]) + kmer['residue'])
            self.calculate_error_weights()

        def solve(self):
            c = config.Configuration()
            if c.solver == 'coin':
                problem, variables = self.generate_mps_linear_program()
                problem.writeLP(os.path.join(self.get_current_job_directory(), 'program_coin_' + str(self.index) + '.lp'))
                problem.writeMPS(os.path.join(self.get_current_job_directory(), 'program_coin_' + str(self.index) + '.mps'))
                command = '/share/hormozdiarilab/Codes/NebulousSerendipity/coin/build/bin/clp ' + os.path.join(self.get_current_job_directory(), 'program_coin_' + str(self.index) + '.mps') + ' -dualsimplex -solution ' + os.path.join(self.get_current_job_directory(), 'solution_' + str(self.index) + '.mps')
                output = subprocess.call(command, shell = True)
                self.import_lp_values('solution_' + str(self.index) + '.mps')
                #with open(os.path.join(self.get_current_job_directory(), 'solution_coin.json'), 'w') as json_file:
                #    json.dump({'variables': self.solution}, json_file, indent = 4, sort_keys = True)
            if c.solver == 'cplex':
                problem = self.generate_linear_program()
                problem.solve()
                problem.write(os.path.join(self.get_current_job_directory(), 'program_cplex.lp'))
                self.solution = problem.solution.get_values()
                with open(os.path.join(self.get_current_job_directory(), 'solution_cplex.json'), 'w') as json_file:
                    json.dump({'variables': self.solution}, json_file, indent = 4, sort_keys = True)
            self.export_solution()

        def generate_mps_linear_program(self):
            c = config.Configuration()
            problem = LpProblem("Nebula", LpMinimize)
            i = 0
            names = [''] * len(self.tracks)
            variables = [None] * (len(self.tracks) + 2 * len(self.lp_kmers))
            regex = re.compile('[^a-zA-Z0-9]')
            for track in self.tracks:
                name = regex.sub('_', track)
                variables[self.tracks[track]['index']] = LpVariable('c' + name, 0, 1)
                problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintLE, 'c' + name + '_ub', 1.0)
                problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintGE, 'c' + name + '_lb', 0.0)
                i += 1
            # error variables
            for index, kmer in enumerate(self.lp_kmers):
                variables[i] = LpVariable('e' + str(index))
                i += 1
            # absolute value of the error variables
            for index, kmer in enumerate(self.lp_kmers):
                variables[i] = LpVariable('l' + str(index))
                i += 1
            expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], kmer['weight']) for index, kmer in enumerate(self.lp_kmers)])
            problem += expr
            for i, kmer in enumerate(self.lp_kmers):
                self.add_mps_error_absolute_value_constraints(problem, variables, i)
                if kmer['type'] == 'junction':
                    indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                    indices.append(len(self.tracks) + i)
                    coeffs = list(map(lambda track: kmer['coverage'][track] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                    coeffs.append(1.0)
                    rhs = kmer['count'] - kmer['residue']
                    expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                    problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
                if kmer['type'] == 'inner':
                    for track in kmer['tracks']:
                        if c.tracks[track].svtype == 'INS':
                            indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                            indices.append(len(self.tracks) + i)
                            coeffs = list(map(lambda track: kmer['coverage'][track] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                            coeffs.append(1.0)
                            rhs = kmer['count'] - kmer['residue']
                            expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                            problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
                        if c.tracks[track].svtype == 'DEL':
                            indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                            indices.append(len(self.tracks) + i)
                            coeffs = list(map(lambda track: -1 * kmer['coverage'][track] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                            coeffs.append(1.0)
                            rhs = kmer['count'] - kmer['residue'] - sum(map(lambda track: kmer['coverage'][track] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                            expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                            problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
                        break
            return problem, variables

        def export_solution(self):
            c = config.Configuration()
            self.errors = self.solution[len(self.tracks):]
            for track in self.tracks:
                index = self.tracks[track]['index']
            self.find_rounding_break_points()
            print('Rounding', len(self.tracks), 'tracks')
            name = 'merge_iteration_' + str(self.index) + '.bed'
            with open(os.path.join(self.get_current_job_directory(), name), 'w') as bed_file:
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
            self.export_kmers()

        def export_kmers(self):
            tracks = {}
            for index, kmer in enumerate(self.lp_kmers):
                for track in kmer['tracks']:
                    if not track in tracks:
                        tracks[track] = {}
                    tracks[track][kmer['kmer']] = kmer
            for track in tracks:
                path = os.path.join(self.get_current_job_directory(), track)
                if not os.path.exists(path):
                    os.makedirs(path)
                with open(os.path.join(path, 'iteration_' + str(self.index) + '.json'), 'w') as json_file:
                    json.dump(tracks[track], json_file, indent = 4)

        def reduce(self):
            c = config.Configuration()
            for track in self.tracks:
                self.tracks[track]['genotypes'] = []
            for iteration in range(self.num_iterations):
                path = os.path.join(self.get_current_job_directory(), 'merge_iteration_' + str(iteration) + '.bed')
                tracks = bed.load_tracks_from_file_as_dict(path)
                for track in tracks:
                    self.tracks[track]['genotypes'].append(tracks[track].lp_genotype)
            with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
                bed_file.write('CHROM\tBEGIN\tEND\tLP_GENOTYPE\tSCORE\tID\n')
                for track in self.tracks:
                    t = c.tracks[track]
                    self.calculate_genotype_score(self.tracks[track])
                    bed_file.write(t.chrom + '\t' +
                        str(t.begin) + '\t' +
                        str(t.end) + '\t' +
                        str(self.tracks[track]['lp_genotype']) + '\t' +
                        str(self.tracks[track]['score']) + '\t' +
                        str(t.id) + '\n')

        def calculate_genotype_score(self, track):
            s = {}
            for g in track['genotypes']:
                if not g in s:
                    s[g] = 0
                s[g] += 1
            m = max(s, key = s.get)
            track['score'] = float(s[m]) / len(track['genotypes'])
            track['lp_genotype'] = m

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExportGenotypingKmersJob(map_reduce.BaseGenotypingJob):

    _name = 'ExportGenotypingKmersJob'
    _category = 'genotyping'
    _previous_job = CgcIntegerProgrammingJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = ExportGenotypingKmersJob(**kwargs)
        job.execute()

    def load_junction_kmer_tracks(self):
        print('selecting junction kmers...')
        self.kmers['junction_kmers'] = {}
        job = CgcJunctionKmersIntegerProgrammingJob()
        tracks = bed.load_tracks_from_file(os.path.join(job.get_current_job_directory(), 'merge.bed'))
        for track in tracks:
            if float(track.lp_value) > 0.10:
                kmers = json.load(open(os.path.join(job.get_current_job_directory(), track.id + '.json'), 'r'))
                self.tracks[track.id] = track
                for kmer in kmers:
                    self.kmers['junction_kmers'][kmer] = kmers[kmer]
                    self.kmers['junction_kmers'][kmer]['count'] = 0
                    self.kmers['junction_kmers'][kmer]['total'] = 0
            else:
                self.not_tracks[track.id] = track

    def load_inner_kmer_tracks(self):
        print('selecting inner kmers...')
        self.kmers['inner_kmers'] = {}
        job = CgcInnerKmersIntegerProgrammingJob()
        tracks = bed.load_tracks_from_file(os.path.join(job.get_current_job_directory(), 'merge.bed'))
        n = []
        m = []
        l = len(self.tracks)
        for track in tracks:
            if float(track.lp_value) > 0.25:
                kmers = json.load(open(os.path.join(job.get_current_job_directory(), track.id + '.json'), 'r'))
                if not track.id in self.tracks: # no junction kmers
                    n.append(track)
                else: # both
                    m.append(track)
                    self.tracks.pop(track.id, None)
                    l -= 1
                for kmer in kmers:
                    self.kmers['inner_kmers'][kmer] = kmers[kmer]
                    self.kmers['inner_kmers'][kmer]['count'] = 0
                    self.kmers['inner_kmers'][kmer]['total'] = 0
                self.not_tracks.pop(track.id, None)
            else:
                self.not_tracks[track.id] = track
        print('selected', l, 'events to genotype with junction kmers only')
        print('selected', len(n), 'events to genotype with inner kmers only')
        print('selected', len(m), 'events to genotype with both inner and junction kmers')
        print('not genotyping', len(self.not_tracks), 'events')
        with open(os.path.join(self.get_current_job_directory(), 'inner.bed'), 'w') as bed_file:
            for track in n:
                bed_file.write(track.serialize())
        with open(os.path.join(self.get_current_job_directory(), 'both.bed'), 'w') as bed_file:
            for track in m:
                bed_file.write(track.serialize())
        with open(os.path.join(self.get_current_job_directory(), 'junction.bed'), 'w') as bed_file:
            for track in self.tracks:
                bed_file.write(self.tracks[track].serialize())

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = {}
        self.not_tracks = {}
        self.load_junction_kmer_tracks()
        self.load_inner_kmer_tracks()
        self._previous_job = CgcCounterJob
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
        self.kmers['depth_kmers'] = kmers['depth_kmers']
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        exit()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcClusteringJob(map_reduce.BaseGenotypingJob):

    _name = 'CgcClusteringJob'
    _category = 'clustering'
    _previous_job = None

    class ClusteringException(Exception):
        pass

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = CgcClusteringJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        #self.paths = ['/share/hormozdiarilab/Codes/NebulousSerendipity/output/clustering/CutlessCandyChr17/' + str(i) + '/CgcIntegerProgrammingJob/merge.bed' for i in range(self.begin, self.end)]
        self.paths = c.bed
        self.load_basis()
        self.tracks = {}
        self.clusters = []
        for index, path in enumerate(self.paths):
            print(path)
            tracks = bed.load_tracks_from_file_as_dict(path, parse_header = True)
            for track in tracks:
                if index == 0 and not track in self.tracks:
                    self.tracks[track] = {}
                if track in self.tracks:
                    self.tracks[track][path] = tracks[track]
        self.round_robin(self.tracks)

    def load_basis(self):
        self.gold = []
        for i, path in enumerate(self.paths):
            sample = path.split('/')[-2]
            self.gold.append(os.path.join('/share/hormozdiarilab/Codes/NebulousSerendipity/data/Audano/Simons', 'Audano.SGDP_' + name + '.hg38.ALL.bed'))

    def transform(self, track, track_name):
        c = config.Configuration()
        features = []
        for index, path in enumerate(self.paths):
            features.append([float(track[path].lp_value)])
        features = np.array(features)
        output_path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        try:
            genotypes, likelihoods = self.fit_gaussian_mixture(track, track_name, features)
            for index, path in enumerate(self.paths):
                track[path].add_field('genotype', genotypes[index])
                track[path].add_field('likelihood', likelihoods[index])
            output_path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        except Exception as e:
            track['error'] = traceback.format_exc()
            print(red(track['error']))
            for index, path in enumerate(self.paths):
                track[path].add_field('genotype', track[path].lp_genotype)
            output_path = os.path.join(self.get_current_job_directory(), track_name + '.error')
            debug_breakpoint()
        with open(output_path, 'w') as json_file:
            json.dump(track, json_file, indent = 4, sort_keys = True, default = bed.BedTrack.json_serialize)
        return output_path

    def select_initial_centroids(self, features, K):
        while True:
            indices = np.random.randint(0, len(features), size = K)
            if len(sorted(set(list(indices)))) == K: #duplicate indices
                tmp = [features[j][0] for j in indices]
                if len(sorted(set(tmp))) == K: #duplicate values
                    break
        centroids = np.array([features[j] for j in indices])
        return centroids

    def find_num_unique_lp_values(self, features):
        n = len(features[0])
        f = [feature[0] for feature in features]
        K = min(3, len(sorted(set(f))))
        return min(3, K)

    def find_num_clusters(self, features):
        #min_lp = np.amin(features)
        #max_lp = np.amax(features)
        f = sorted([feature[0] for feature in features])
        p = f[0]
        k = 1
        for c in f[1:]:
            if c - p > 0.15:
                k += 1
            p = c
        self.clusters.append(k)
        return self.find_num_unique_lp_values(features)
        if k >= 3:
            return 3
        min_lp = np.amin(features)
        max_lp = np.amax(features)
        median_lp = np.median(features)
        k = 0
        if min_lp <= 0.1:
            k += 1
        if max_lp >= 0.8:
            k += 1

    def fit_gaussian_mixture(self, track, track_name, features):
        print(cyan(track_name))
        means = [0.0, 0.0, 0.0]
        max_likelihood = -1
        choice = None
        k = self.find_num_unique_lp_values(features)
        if k == 1:
            m = np.argmax(features)
            print(yellow('HERE'))
            return [track[self.paths[m]].lp_genotype] * len(features), [1.0] * len(features)
        features = features.flatten()
        for K in range(2, k + 1):
            mean_step = 0.05
            std_step = 0.03
            if K == 1:
                mean_choices = [(mean_step * i,) for i in range(0, 3)]
                std_choices = [(std_step * i,) for i in range(1, 6)]
            if K == 2:
                mean_choices = [(mean_step * i, mean_step * j) for i in range(0, int(1.0 / mean_step) - 4) for j in range(i + 3, int(1.0 / mean_step) + 1)]
                std_choices = [(std_step * i, std_step * j) for i in range(1, 6) for j in range(1, 5)]
            if K == 3:
                mean_choices = [(mean_step * 0, mean_step * j, mean_step * k) for j in range(4, int(1.0 / mean_step) - 4) for k in range(j + 3, int(1.0 / mean_step) + 1)]
                std_choices = [(std_step * i, std_step * j, std_step * k) for i in range(1, 6) for j in range(1, 6) for k in range(1, 5)]
            for mean_choice in mean_choices:
                for std_choice in std_choices:
                    #print(mean_choice)
                    #print(std_choice)
                    l = self.calculate_likelihood(features, mean_choice, std_choice)
                    if l > max_likelihood:
                        max_likelihood = l
                        choice = (mean_choice, std_choice)
        print(yellow(choice))
        self.plot_distribution(features, choice[0], choice[1], track_name)
        #debug_breakpoint()
        return self.likelihood_genotype(features, choice)

    #def fit_normal_dsitributions_em(self, trakc, track_name, features):
    #    for K in range(2, 3 + 1):

    def likelihood_genotype(self, features, choice):
        distributions = [statistics.NormalDistribution(mean, std) for mean, std in zip(choice[0], choice[1])]
        genotypes = ['00', '10', '11']
        likelihoods = [np.argmax([distribution.pmf(feature) for distribution in distributions]) for feature in features]
        return [genotypes[i] for i in likelihoods], [distributions[i].pmf(feature) for feature, i in zip(features, likelihoods)]

    def calculate_likelihood(self, features, means, stds):
        distributions = [statistics.NormalDistribution(mean, std) for mean, std in zip(means, stds)]
        likelihood = sum([max([distribution.pmf(feature) for distribution in distributions]) for feature in features])
        return likelihood

    def plot_distribution(self, features, mean, std, track_name):
        #x = [0.05 * i for i in range(0, 21)]
        #data = [graph_objs.Scatter(x = x, y = stats.norm(loc = m, scale = s).pdf(x) * 10, mode = 'lines', line = dict(color = 'rgba(0, 0, 0)')) for m, s in zip(mean, std)]
        #data.append(graph_objs.Histogram(x = features, xbins = dict(start = 0.0, size = 0.05, end = 1.0)))
        #layout = graph_objs.Layout(title = track_name)
        #figure = graph_objs.Figure(data = data, layout = layout)
        #plotly.plot(figure, filename = os.path.join(self.get_current_job_directory(), 'normal_overlay_' + track_name + '.html'), auto_open = False)
        pass

    def kmeans(self, track, track_name, features):
        print(blue('clustering', track_name))
        K = self.find_num_clusters(features)
        if K == 1:
            m = np.argmax(features)
            return [track[self.paths[m]].lp_genotype] * len(features), [1.0] * len(features)
        centroids = self.select_initial_centroids(features, K)
        old_centroids = np.zeros(centroids.shape)
        candidates = []
        errors = []
        error = 0
        m = 0
        r = 0
        n = 0
        #print('round 0 centroids:', centroids)
        while True:
            m += 1
            n += 1
            if m == 100 or m == -1 or (old_centroids == centroids).all():
                if m == -1:
                    print(red('empty cluster, skipping round'))
                if m == 100:
                    print(red('didn\'t converge, skipping round'))
                    debug_breakpoint()
                else:
                    #print(green('round', r), cyan('iteration', m))
                    m = 0
                    errors.append(error)
                    candidates.append(centroids)
                    centroids = self.select_initial_centroids(features, K)
                    #print('round', r, 'centroids:', centroids)
                    old_centroids = np.zeros(centroids.shape)
                r += 1
                if r == 10:
                    break
            clusters, error = self.cluster_data(features, centroids)
            #print('updated centroids:', centroids)
            old_centroids = copy.deepcopy(centroids)
            for i in range(K):
                points = np.array([features[j] for j in range(len(features)) if clusters[j] == i])
                #print('cluster', i, 'points', points)
                if len(points) == 0:
                    m = -1
                    print(red('invalid clustering'))
                    continue
                centroids[i] = np.mean(points, axis = 0)
        choice = np.argmin(errors)
        centroids = candidates[choice]
        clusters, _ = self.cluster_data(features, centroids)
        genotypes = self.assign_genotypes(features, clusters, centroids)
        self.plot_clusters(track_name, features, clusters)
        return genotypes, [1.0] * len(features)

    def cluster_data(self, features, centroids):
        error = 0
        clusters = []
        for i in range(len(features)):
            distances = self.distance(features[i], centroids)
            cluster = np.argmin(distances)
            error += distances[cluster]
            clusters.append(cluster)
        return clusters, error

    def assign_genotypes(self, features, clusters, centroids):
        K = len(centroids)
        min_lp = np.amin(features)
        max_lp = np.amax(features)
        g = []
        if K == 3:
            g = ['00', '10', '11']
        if K == 2:
            if max_lp > 0.75:
                if min_lp > 0.1:
                    g = ['10', '11']
                else:
                    g = ['00', '11']
            else:
                g = ['00', '10']
        if K == 1:
            if max_lp > 0.75:
                g = ['11']
            elif max_lp > 0.25:
                g = ['10']
            else:
                g = ['00']
        # this part is highly heuristic
        #print(centroids)
        genotypes = []
        norms = np.array([np.linalg.norm(centroid) for centroid in centroids])
        centroids = []
        for i in range(len(norms)):
            centroid = np.argmin(norms)
            norms[centroid] += 2 # To make it bigger than whatever other centroid remaining
            centroids.append(centroid)
        #print(centroids)
        #print(g)
        #print(zip(clusters, features))
        for i in range(len(features)):
            for j in range(len(centroids)):
                if clusters[i] == centroids[j]:
                    genotypes.append(g[j])
        #print(zip(clusters, genotypes))
        return genotypes

    def plot_clusters(self, track, features, clusters):
        x = []
        y = []
        for i in range(len(features)):
            x.append(clusters[i])
            y.append(features[i][0])
        try:
            visualizer.histogram(y, 'kmeans_' + track, self.get_current_job_directory(), 'cluster', 'value', step = 0.05)
            visualizer.violin(x, y, 'kmeans_' + track, self.get_current_job_directory(), 'cluster', 'value')
        except Exception as e:
            print(yellow(e))

    def distance(self, a, b):
        return np.linalg.norm(a - b, axis = None if a.shape == b.shape else 1)

    def reduce(self):
        print('reducing')
        files = {}
        for path in self.paths:
            name = path.split('/')[-1]
            files[path] = open(os.path.join(self.get_current_job_directory(), name), 'w')
            files[path].write('CHROM\tBEGIN\tEND\tLP_GENOTYPE\tLP_VALUE\tLIKELIHOOD\tID\n')
        for batch in self.load_output():
            for track in batch:
                track = json.load(open(batch[track]))
                for path in track:
                    if path == 'error':
                        continue
                    t = track[path]
                    files[path].write(t['chrom'] + '\t' + str(t['begin']) + '\t' + str(t['end']) + '\t' + t['genotype'] + '\t' + t['lp_value'] + '\t' + str(t['likelihood']) + '\t' + t['id'] + '\n')
        #self.gather_genotype_statistics()

    def simulation_analysis(self):
        files = {}
        for path in self.paths:
            name = path.split('/')[-1]
            files[path] = open(os.path.join(os.path.split(path)[0], 'cluster.bed'), 'w')
            files[path].write('CHROM\tBEGIN\tEND\tLP_GENOTYPE\tLP_VALUE\tLIKELIHOOD\tID\n')
        for batch in self.load_output():
            for track in batch:
                track = json.load(open(batch[track]))
                for path in track:
                    if path == 'error':
                        continue
                    t = track[path]
                    files[path].write(t['chrom'] + '\t' + str(t['begin']) + '\t' + str(t['end']) + '\t' + t['genotype'] + '\t' + t['lp_value'] + '\t' + str(t['likelihood']) + '\t' + t['id'] + '\n')
        self.tracks = {}

    def gather_genotype_statistics(self):
        c = config.Configuration()
        self.stats = {}
        self.tracks = {}
        for name in ['merge', 'cluster']:
            self.stats[name] = {'genotype': {}, 'likelihood': {}, 'state': {}}
        for p in ['00', '10', '11']:
            for q in ['00', '10', '11']:
                for name in ['merge', 'cluster']:
                    s = p + '_as_' + q
                    self.stats[name]['genotype'][s] = 0
        for p in ['t', 'f']:
            for q in ['p', 'n']:
                for name in ['merge', 'cluster']:
                    self.stats[name]['state'][p + q] = 0
                    self.stats[name]['likelihood'][p + q] = []
        cwd = c.workdir
        for i in range(len(self.paths)):
        #for i in range(self.begin, self.end):
            #cwd = '/share/hormozdiarilab/Codes/NebulousSerendipity/output/clustering/CutlessCandyChr17/' + str(i) + '/CgcIntegerProgrammingJob'
            if i == 1000:
                tracks = bed.load_tracks_from_file_as_dict(self.paths[i])
                #tracks = bed.load_tracks_from_file_as_dict(os.path.join(cwd, 'merge.bed'))
                for track in tracks:
                    self.tracks[track] = {'merge': {}, 'cluster': {}}
                for p in ['00', '10', '11']:
                    for q in ['00', '10', '11']:
                        for track in self.tracks:
                            self.tracks[track]['merge'][p + '_as_' + q] = [0, []]
                            self.tracks[track]['cluster'][p + '_as_' + q] = [0, []]
            self.tabulate('merge', cwd, i)
            self.tabulate('cluster', cwd, i)
        self.plot_likelihoods()
        with open(os.path.join(self.get_current_job_directory(), 'statistics.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'comparison.json'), 'w') as json_file:
            json.dump(self.stats, json_file, indent = 4)

    def plot_likelihoods(self):
        #data = [graph_objs.Histogram(x = self.stats['cluster']['likelihood'][s]) for s in ['tp', 'tn', 'fp', 'fn']]
        #layout = graph_objs.Layout(title = 'Likelihoods')
        #figure = graph_objs.Figure(data = data, layout = layout)
        #plotly.plot(figure, filename = os.path.join(self.get_current_job_directory(), 'likelihoods.html'), auto_open = False)
        pass

    def tabulate(self, name, cwd, i):
        p = subprocess.Popen(['/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/tabulate.sh', name + '.bed'], cwd = cwd)
        p.wait()
        #p = subprocess.Popen(['/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/verify_sim.sh', '/share/hormozdiarilab/Codes/NebulousSerendipity/output/simulation/CutlessCandyChr17/{}/Simulation'.format(str(i))], cwd = cwd)
        p = subprocess.Popen(['/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/verify_sim.sh', self.gold[i]], cwd = cwd)
        p.wait()
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        for p in ['00', '10', '11']:
            for q in ['00', '10', '11']:
                s = p + '_as_' + q
                tracks = bed.load_tracks_from_file_as_dict(os.path.join(cwd, s + '.bed'))
                for track in tracks:
                    self.tracks[track][name][s][0] += 1
                    self.tracks[track][name][s][1].append(i)
                    self.stats[name]['genotype'][s] += 1
                    if p == q:
                        if '1' in p:
                            tp += 1
                            self.stats[name]['state']['tp'] += 1
                            if name == 'cluster':
                                self.stats[name]['likelihood']['tp'].append(tracks[track]['likelihood'])
                        else:
                            tn += 1
                            self.stats[name]['state']['tn'] += 1
                            if name == 'cluster':
                                self.stats[name]['likelihood']['tn'].append(tracks[track]['likelihood'])
                    else:
                        if not '1' in p:
                            fp += 1
                            self.stats[name]['state']['fp'] += 1
                            if name == 'cluster':
                                self.stats[name]['likelihood']['fp'].append(tracks[track]['likelihood'])
                        else:
                            if '1' in q:
                                tp += 1
                                self.stats[name]['state']['tp'] += 1
                                if name == 'cluster':
                                    self.stats[name]['likelihood']['tp'].append(tracks[track]['likelihood'])
                            else:
                                fn += 1
                                self.stats[name]['state']['fn'] += 1
                                if name == 'cluster':
                                    self.stats[name]['likelihood']['fn'].append(tracks[track]['likelihood'])

    def kmer_distribution(self):
        kmer = 'CCCCGCTATTATTTCCCAGGTAGCTGGGACTA' # HG00514_chr17-617898-INS-50
        kmer = 'GATCTAATTAAATTAAAGAGCTTCTGTACAGC' # CHM1_chr17-6944085-DEL-736
        kmer = 'CTCTCCCTCTCCTCTCCCTCTCTCCCTCTCCC' # HG00514_chr17-20038316-INS-2590
        kmer = 'CGGGCGGCTGGCCGGGCGGGGGCTGACCCCCA'
        counts = []
        for i in range(self.begin, self.end):
            cwd = '/share/hormozdiarilab/Codes/NebulousSerendipity/output/clustering/CutlessCandyChr17/' + str(i) + '/CgcCounterJob'
            print(cwd)
            with open(os.path.join(cwd, 'HG00514_chr17-20038316-INS-2590.json'), 'r') as json_file:
                kmers = json.load(json_file)
                counts.append(kmers['inner_kmers'][kmer]['count'] if kmer in kmers['inner_kmers'] else kmers['junction_kmers'][kmer]['count'])
        visualizer.histogram(counts, kmer, self.get_current_job_directory(), 'count', 'samples')
