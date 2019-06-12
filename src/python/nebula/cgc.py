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
    gapped,
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

import plotly.offline as plotly
import plotly.graph_objs as graph_objs

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CgcCounterJob(map_reduce.FirstGenotypingJob, counter.BaseExactCountingJob):

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
        self.gapped_kmers = {}
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
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            kmers = {}
            kmers['junction_kmers'] = self.junction_kmers
            kmers['inner_kmers'] = self.inner_kmers
            kmers['depth_kmers'] = self.depth_kmers
            kmers['gc_kmers'] = self.gc_kmers
            json.dump(kmers, json_file, indent = 4)

    def export_counter_input(self):
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
        if kmer in self.gapped_kmers:
            count = tokens[0] 
            self.gapped_kmers[kmer]['count'] += count
        else:
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
        #counts = [[self.gc_kmers[kmer]['count'] for kmer in list(filter(lambda k: self.gc_kmers[kmer]['gc'] == gc, self.gc_kmers))], range(0, 96 + 1)]
        #mean = np.mean(counts, axis = 1)
        #std = np.std(counts, axis = 1)

    def reduce(self):
        c = config.Configuration()
        self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'depth_kmers'), 'w') as json_file:
            json.dump(self.depth_kmers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
            json.dump(self.inner_kmers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
            json.dump(self.gapped_kmers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'junction_kmers.json'), 'w') as json_file:
            json.dump(self.junction_kmers, json_file, indent = 4)
        if c.cgc:
            kmers = {}
            kmers['depth_kmers'] = self.depth_kmers
            kmers['inner_kmers'] = self.inner_kmers
            kmers['junction_kmers'] = self.junction_kmers
            name = 'kmers_' + (c.fastq.split('/')[-1] if c.fastq else c.bam.split('/')[-1]) + '.json'
            with open(os.path.join(os.getcwd(), name), 'w') as json_file:
                json.dump(kmer, json_file, indent = 4)
        self.tracks = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.gapped_kmers:
            for track in self.gapped_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['gapped_kmers'][kmer] = self.gapped_kmers[kmer]
                self.tracks[track]['gapped_kmers'][kmer]['reference'] = 1
                self.tracks[track]['gapped_kmers'][kmer]['type'] = 'gapped'
        for kmer in self.junction_kmers:
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)
        return self.estimate_depth_of_coverage()

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
        self.round_robin(self.load_previous_job_results())
        self.resume_from_reduce = False
        self.lp_kmers = {}

    def transform(self, path, track_name):
        print(green(track_name))
        kmers = json.load(open(os.path.join(self.get_previous_job_directory(), path), 'r'))
        c = config.Configuration()
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
        for kmer in kmers['gapped_kmers']:
            lp_kmers[kmer] = True
            self.lp_kmers[kmer] = kmers['gapped_kmers'][kmer]
            self.lp_kmers[kmer]['reference'] = 1 
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        print(path)
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4, sort_keys = True)
        return path

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
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
        gapped = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'gapped']
        junction = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'junction']
        for kmers in [inner, gapped, junction]:
            for kmer in kmers:
                kmer['weight'] = float(len(inner) + len(gapped) + len(junction)) / (3.0 * len(kmers))

    def generate_linear_program(self):
        c = config.Configuration()
        #globals()['cplex'] = __import__('cplex')
        #problem = cplex.Cplex()
        #problem.objective.set_sense(problem.objective.sense.minimize)
        ## the coverage of each event
        #names = [''] * len(self.tracks)
        #regex = re.compile('[^a-zA-Z0-9]')
        #for track in self.tracks:
        #    name = regex.sub('', track)
        #    tokens = name.split('_')
        #    names[self.tracks[track]['index']] = 'c' + name
        #problem.variables.add(names = names,
        #    ub = [1.0] * len(self.tracks),
        #)
        ## the real-valued error parameter for kmer
        #problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
        #    lb = [-100000000 for kmer in self.lp_kmers]
        #)
        ## absolute value of the kmer error parameter
        #problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
        #    obj = [kmer['weight'] for index, kmer in enumerate(self.lp_kmers)]
        #)
        #print('adding constraints')
        #for index, kmer in enumerate(self.lp_kmers):
        #    self.add_error_absolute_value_constraints(problem, index)
        #    if kmer['type'] == 'junction':
        #        ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
        #        ind.append(len(self.tracks) + index) # Objective
        #        val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
        #        val.append(1.0) #Objective
        #        rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']]
        #        problem.linear_constraints.add(
        #            lin_expr = [cplex.SparsePair(
        #                ind = ind,
        #                val = val,
        #            )],
        #            rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']],
        #            senses = ['E']
        #        )
        #    if kmer['type'] == 'inner':
        #        for track in kmer['tracks']:
        #            if c.tracks[track].svtype == 'INS':
        #                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
        #                ind.append(len(self.tracks) + index) # Objective
        #                val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
        #                val.append(1.0) #Objective
        #                rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']]
        #                problem.linear_constraints.add(
        #                    lin_expr = [cplex.SparsePair(
        #                        ind = ind,
        #                        val = val,
        #                    )],
        #                    rhs = rhs,
        #                    senses = ['E']
        #                )
        #            else:
        #                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
        #                ind.append(len(self.tracks) + index) # Objective
        #                val = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
        #                val.append(1.0) #Objective
        #                rhs = [kmer['count'] - kmer['coverage'] * kmer['residue'] - sum(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))]
        #                problem.linear_constraints.add(
        #                    lin_expr = [cplex.SparsePair(
        #                        ind = ind,
        #                        val = val,
        #                    )],
        #                    rhs = rhs,
        #                    senses = ['E']
        #                )
        #            break
        #print('completed program generation')
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
            if kmer['type'] == 'junction':
                indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                indices.append(len(self.tracks) + i)
                coeffs = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['count'] - kmer['coverage'] * kmer['residue']
                expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
            if kmer['type'] == 'inner':
                for track in kmer['tracks']:
                    if c.tracks[track].svtype == 'INS':
                        indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                        indices.append(len(self.tracks) + i)
                        coeffs = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                        coeffs.append(1.0)
                        rhs = kmer['count'] - kmer['coverage'] * kmer['residue']
                        expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                        problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
                    if c.tracks[track].svtype == 'DEL':
                        indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                        indices.append(len(self.tracks) + i)
                        coeffs = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                        coeffs.append(1.0)
                        rhs = kmer['count'] - kmer['coverage'] * kmer['residue'] - sum(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                        expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                        problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
                    break
        return problem, variables

    #def create_output_directories(self):
    #    c = config.Configuration()
    #    if not c.cgc:
    #        map_reduce.Job.create_output_directories(self)

    #def get_current_job_directory(self):
    #    c = config.Configuration()
    #    if c.cgc:
    #        return os.getcwd()
    #    return map_reduce.BaseGenotypingJob.get_current_job_directory(self)

    def round_genotype(self, c, svtype):
        if c > self.b2:
            return (1.0, '11')
        elif c > self.b1:
            return (0.5, '10')
        else:
            return (0.0, '00')

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

class CgcCoverageCorrectingIntegerProgrammingJob(CgcInnerKmersIntegerProgrammingJob):

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
            self.coverages = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
            kmer['coverage'] = {track: np.random.normal(c.coverage, c.std) for track in kmer['tracks']}
            #kmer['coverage'] = {track: self.coverages[random.randint(0, len(self.coverages) - 1)] for track in kmer['tracks']}
            kmer['residue'] = sum([np.random.normal(c.coverage, c.std) for i in range(kmer['reference'] - r)])
            #kmer['residue'] = sum([self.coverages[random.randint(0, len(self.coverages) - 1)] for i in range(kmer['reference'] - r)])
            kmer['weight'] = 1.0
            kmer['count'] = min(kmer['count'], sum(kmer['coverage'].values()) + kmer['residue'])
        self.calculate_error_weights()

    def solve(self):
        c = config.Configuration()
        self.num_iterations = 10
        for iteration in range(self.num_iterations):
            self.calculate_residual_coverage()
            if c.solver == 'coin':
                problem, variables = self.generate_mps_linear_program()
                problem.writeLP(os.path.join(self.get_current_job_directory(), 'program_coin_' + str(iteration) + '.lp'))
                problem.writeMPS(os.path.join(self.get_current_job_directory(), 'program_coin_' + str(iteration) + '.mps'))
                command = '/share/hormozdiarilab/Codes/NebulousSerendipity/coin/build/bin/clp ' + os.path.join(self.get_current_job_directory(), 'program_coin_' + str(iteration) + '.mps') + ' -dualsimplex -solution ' + os.path.join(self.get_current_job_directory(), 'solution.mps')
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
            self.export_solution(iteration)
        return self.tracks

    def export_solution(self, iteration):
        c = config.Configuration()
        self.errors = self.solution[len(self.tracks):]
        for track in self.tracks:
            index = self.tracks[track]['index']
            self.tracks[track]['coverage'].append(self.solution[index])
            self.tracks[track]['lp_kmers'] = []
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['lp_kmers'].append(kmer)
        self.find_rounding_break_points()
        print('Rounding', len(self.tracks), 'tracks')
        name = 'merge_iteration_' + str(iteration) + '.bed'
        with open(os.path.join(self.get_current_job_directory(), name), 'w') as bed_file:
            bed_file.write('CHROM\tBEGIN\tEND\tLP_GENOTYPE\tLP_VALUE\tSCORE\tID\n')
            for track in self.tracks:
                t = c.tracks[track]
                index = self.tracks[track]['index']
                g = self.round_genotype(self.solution[index], t.svtype)
                self.tracks[track]['genotypes'].append(g[1])
                self.calculate_genotype_score(self.tracks[track])
                bed_file.write(t.chrom + '\t' +
                    str(t.begin) + '\t' +
                    str(t.end) + '\t' +
                    str(self.tracks[track]['lp_genotype']) + '\t' +
                    str(self.solution[index]) + '\t' + 
                    str(self.tracks[track]['score']) + '\t' +
                    str(t.id) + '\n')
        self.export_kmers_for_iteration(iteration)
        if iteration == self.num_iterations - 1:
            for track in self.tracks:
                visualizer.violin([1] * self.num_iterations, self.tracks[track]['coverage'], track, self.get_current_job_directory(), track, 'LP')

    def export_kmers_for_iteration(self, iteration):
        tracks = {}
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                if not track in tracks:
                    tracks[track] = {}
                tracks[track][kmer['kmer']] = kmer
        for track in tracks:
            path = os.path.join(self.get_current_job_directory(), track)
            if iteration == 0:
                if not os.path.exists(path):
                    os.makedirs(path)
            with open(os.path.join(path, 'iteration_' + str(iteration) + '.json'), 'w') as json_file:
                json.dump(tracks[track], json_file, indent = 4)

    def calculate_genotype_score(self, track):
        s = {}
        for g in track['genotypes']:
            if not g in s:
                s[g] = 0
            s[g] += 1
        m = max(s, key = s.get)
        track['score'] = float(s[m]) / len(track['genotypes'])
        track['lp_genotype'] = m

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
        self.kmers['half_mers'] = {}
        self.kmers['depth_kmers'] = kmers['depth_kmers']
        self.kmers['gapped_kmers'] = {}
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        exit()

    #def load_inputs(self):
    #    c = config.Configuration()
    #    self.tracks = {}
    #    tracks = bed.load_tracks_from_file_as_dict(os.path.join(self.get_previous_job_directory(), 'merge.bed'))
    #    for track in tracks:
    #        if float(tracks[track].lp_value) < 0.05:
    #            print(track)
    #            self.tracks[track] = True
    #    self._previous_job = CgcCounterJob
    #    with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
    #        self.kmers = json.load(json_file)
    #    self.half_mers = {kmer: self.kmers['half_mers'][kmer] for kmer in self.kmers['half_mers'] if not self.has_track_overlap(self.kmers['half_mers'][kmer])}
    #    self.depth_kmers = self.kmers['depth_kmers']
    #    self.inner_kmers = {kmer: self.kmers['inner_kmers'][kmer] for kmer in self.kmers['inner_kmers'] if not self.has_track_overlap(self.kmers['inner_kmers'][kmer])}
    #    self.gapped_kmers = {kmer: self.kmers['gapped_kmers'][kmer] for kmer in self.kmers['gapped_kmers'] if not self.has_track_overlap(self.kmers['gapped_kmers'][kmer])}
    #    print(len(self.kmers['junction_kmers']))
    #    self.junction_kmers = {kmer: self.kmers['junction_kmers'][kmer] for kmer in self.kmers['junction_kmers'] if not self.has_track_overlap(self.kmers['junction_kmers'][kmer])}
    #    print(len(self.junction_kmers))
    #    with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
    #        kmers = {}
    #        kmers['junction_kmers'] = self.junction_kmers
    #        kmers['gapped_kmers'] = self.gapped_kmers
    #        kmers['inner_kmers'] = self.inner_kmers
    #        kmers['depth_kmers'] = self.depth_kmers
    #        kmers['half_mers'] = self.half_mers
    #        json.dump(kmers, json_file, indent = 4)
    #    exit()

    def has_track_overlap(self, a):
        for t in a['tracks']:
            if t in self.tracks:
                return True
        return False

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
        self.paths = c.bed
        self.tracks = {}
        self.clusters = []
        for index, path in enumerate(c.bed):
            print(path)
            tracks = bed.load_tracks_from_file_as_dict(path, parse_header = True)
            for track in tracks:
                if index == 0 and not track in self.tracks:
                    self.tracks[track] = {}
                if track in self.tracks:
                    self.tracks[track][path] = tracks[track]
        self.round_robin(self.tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        features = []
        for index, path in enumerate(self.paths):
            features.append([float(track[path].lp_value)])
        features = np.array(features)
        try:
            genotypes = self.kmeans(track, track_name, features)
            for index, path in enumerate(self.paths):
                track[path].add_field('genotype', genotypes[index])
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
        return min(3, k)

    def kmeans(self, track, track_name, features):
        print(blue('clustering', track_name))
        K = self.find_num_clusters(features)
        if K == 1:
            m = np.argmax(features)
            return [track[self.paths[m]].lp_genotype] * len(features)
        centroids = self.select_initial_centroids(features, K)
        old_centroids = np.zeros(centroids.shape)
        candidates = []
        errors = []
        error = 0
        m = 0
        r = 0
        print('round 0 centroids:', centroids)
        while True:
            m += 1
            if m == 100 or m == -1 or (old_centroids == centroids).all():
                if m == -1:
                    print(red('empty cluster, skipping round'))
                if m == 100:
                    print(red('didn\'t converge, skipping round'))
                    debug_breakpoint()
                else:
                    print(green('round', r), cyan('iteration', m))
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
        return genotypes

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
        # this part is highly heuristic
        else:
            if min_lp > 0.2:
                g = ['10', '11']
            else:
                g = ['00', '10']
        genotypes = []
        norms = np.array([np.linalg.norm(centroid) for centroid in centroids])
        centroids = []
        for i in range(len(norms)):
            centroid = np.argmin(norms)
            norms[centroid] += 1
            centroids.append(centroid)
        #if len(g) != len(centroids):
        #    print(centroids)
        #    print(g)
        #    print(clusters)
        for i in range(len(features)):
            for j in range(len(centroids)):
                if clusters[i] == centroids[j]:
                    genotypes.append(g[j])
        return genotypes

    #def find_minimum_nudge(self, features, centroids, clusters):
    #    m = 10000
    #    for i in range(len(features)):
    #        c = centroids[clusters[i]]
    #        d = self.distance(features[i], c)
    #        n = [self.distance(features[i], centroid) - d for centroid in centroids if centroid != c]
    #        m = min(m, min(n))
    #    return m

    def plot_clusters(self, track, features, clusters):
        x = []
        y = []
        for i in range(len(features)):
            x.append(clusters[i])
            y.append(features[i][0])
        try:
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
        for batch in self.load_output():
            for track in batch:
                print(batch[track])
                track = json.load(open(batch[track]))
                for path in track:
                    if path == 'error':
                        continue
                    t = track[path]
                    files[path].write(t['chrom'] + '\t' + str(t['begin']) + '\t' + str(t['end']) + '\t' + t['genotype'] + '\t' + t['id'] + '\t' + t['lp_value'] + '\n')

