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
    #visualizer,
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
    _category = 'preprocessing'
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
        self.half_mers = self.kmers['half_mers']
        self.depth_kmers = self.kmers['depth_kmers']
        self.inner_kmers = self.kmers['inner_kmers']
        self.gapped_kmers = self.kmers['gapped_kmers']
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
        print('Gapped kmers:', blue(len(self.gapped_kmers)))
        print('Depth kmers:', blue(len(self.depth_kmers)))
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            kmers = {}
            kmers['junction_kmers'] = self.junction_kmers
            kmers['gapped_kmers'] = self.gapped_kmers
            kmers['inner_kmers'] = self.inner_kmers
            kmers['depth_kmers'] = self.depth_kmers
            kmers['half_mers'] = self.half_mers
            json.dump(kmers, json_file, indent = 4)

    def export_counter_input(self):
        with open(os.path.join(self.get_current_job_directory(), 'pre_inner_kmers.json'), 'w') as json_file:
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
            print('Dumping', green(len(_kmers)), 'kmers for the counter')
            json.dump(_kmers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
            json.dump(self.half_mers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'pre_gapped_kmers.json'), 'w') as json_file:
            json.dump(self.gapped_kmers, json_file, indent = 4)

    def merge_count(self, kmer, tokens):
        if kmer in self.gapped_kmers:
            count = tokens[0] 
            self.gapped_kmers[kmer]['count'] += count
        else:
            count = tokens[0] 
            total = tokens[1] 
            canon = canonicalize(kmer)
            if canon in self.inner_kmers:
                self.inner_kmers[canon]['count'] += count / 2
                self.inner_kmers[canon]['total'] += total / 2
            elif canon in self.depth_kmers:
                self.depth_kmers[canon]['count'] += count / 2
            else:
                self.junction_kmers[canon]['count'] += count / 2
                self.junction_kmers[canon]['total'] += total / 2

    def estimate_depth_of_coverage(self):
        c = config.Configuration()
        self.counts = list(map(lambda kmer: self.depth_kmers[kmer]['count'], self.depth_kmers))
        self.mean = np.mean(self.counts)
        self.std = np.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter outliers
        self.counts = list(filter(lambda x: x < 3 * self.mean, self.counts))
        self.mean = np.mean(self.counts)
        self.std = np.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter outliers
        self.counts = list(filter(lambda x: x < 2 * self.mean, self.counts))
        self.mean = np.mean(self.counts)
        self.std = np.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        #
        stats = {'coverage': self.mean, 'std': self.std}
        with open(os.path.join(self.get_current_job_directory(), 'stats_' + str(c.ksize) + '.json'), 'w') as json_file:
            json.dump(stats, json_file, sort_keys = True, indent = 4)
        return stats

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
    _category = 'preprocessing'
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
            #l = len(self.lp_kmers)
            #l = len(list(filter(lambda i: self.lp_kmers[i]['type'] != 'gapped', self.tracks[kmer['tracks'].keys()[0]]['kmers'])))
            #l = l if l else 1
            #kmer['weight'] = l if kmer['type'] == 'gapped' else 1.0
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
        pass
        #c = config.Configuration()
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
        #return problem

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
                    else:
                        indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                        indices.append(len(self.tracks) + i)
                        coeffs = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                        coeffs.append(1.0)
                        rhs = kmer['count'] - kmer['coverage'] * kmer['residue'] - sum(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                        expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                        problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
                    break
        return problem, variables

    def create_output_directories(self):
        c = config.Configuration()
        if not c.cgc:
            map_reduce.Job.create_output_directories(self)

    def get_current_job_directory(self):
        c = config.Configuration()
        if c.cgc:
            return os.getcwd()
        return map_reduce.BaseGenotypingJob.get_current_job_directory(self)

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
    _category = 'preprocessing'
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
            if track.lp_value > 0.10:
                kmers = json.load(open(os.path.join(job.get_current_job_directory(), track.id + '.json'), 'r'))
                self.tracks[track.id] = kmers
                for kmer in kmers:
                    self.kmers['junction_kmers'][kmer] = kmers[kmer]

    def load_inner_kmer_tracks(self):
        print('selecting inner kmers...')
        self.kmers['inner_kmers'] = {}
        job = CgcInnerKmersIntegerProgrammingJob()
        tracks = bed.load_tracks_from_file(os.path.join(job.get_current_job_directory(), 'merge.bed'))
        for track in tracks:
            if track.lp_value > 0.25:
                kmers = json.load(open(os.path.join(job.get_current_job_directory(), track.id + '.json'), 'r'))
                if not track.id in self.tracks: # no junction kmers
                    self.tracks[track.id] = kmers 
                else: # both
                    self.tracks[track.id].update(kmers)
                for kmer in kmers:
                    self.kmers['inner_kmers'][kmer] = kmers[kmer]

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = {}
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

#class CgcClusteringJob(map_reduce.BaseGenotypingJob):
#
#    _name = 'CgcClusteringJob'
#    _category = 'clustering'
#    _previous_job = None
#
#    class ClusteringException(Exception):
#        pass
#
#    # ============================================================================================================================ #
#    # Launcher
#    # ============================================================================================================================ #
#
#    @staticmethod
#    def launch(**kwargs):
#        c = config.Configuration()
#        job = CgcClusteringJob(**kwargs)
#        job.execute()
#
#    # ============================================================================================================================ #
#    # MapReduce overrides
#    # ============================================================================================================================ #
#
#    def load_inputs(self):
#        c = config.Configuration()
#        self.paths = c.bed
#        self.tracks = {}
#        for index, path in enumerate(c.bed):
#            tracks = bed.load_tracks_from_file_as_dict(path, parse_header = True)
#            for track in tracks:
#                if index == 0 and not track in self.tracks:
#                    self.tracks[track] = {}
#                if track in self.tracks:
#                    self.tracks[track][path] = tracks[track]
#        self.round_robin(self.tracks)
#
#    def transform(self, track, track_name):
#        c = config.Configuration()
#        features = []
#        for index, path in enumerate(self.paths):
#            features.append([float(track[path].lp_value)])
#        features = np.array(features)
#        genotypes, centroids = self.kmeans(track_name, features)
#        for index, path in enumerate(self.paths):
#            track[path].add_field('genotype', genotypes[index])
#        t = {}
#        for index, path in enumerate(self.paths):
#            t[path] = track[path].as_dict()
#        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
#        with open(path, 'w') as json_file:
#            json.dump(t, json_file, indent = 4, sort_keys = True)
#        return path
#
#    def find_num_clusters(self, features):
#        n = len(features[0])
#        f = [feature[0] for feature in features]
#        K = min(3, len(sorted(set(f))))
#        if K != 3:
#            raise CgcClusteringJob.ClusteringException
#        return K
#
#    def select_initial_centroids(self, features, K):
#        while True:
#            indices = np.random.randint(0, len(features), size = K)
#            print(indices)
#            if len(sorted(set(list(indices)))) == K: #duplicate indices
#                tmp = [features[j][0] for j in indices]
#                if len(sorted(set(tmp))) == K: #duplicate values
#                    break
#        centroids = np.array([features[j] for j in indices])
#        return centroids
#
#    def kmeans(self, track, features):
#        print(blue('clustering', track))
#        #if track != 'CHM1_chr1-100528655-INS-331':
#        #    raise CgcCounterJob.ClusteringException
#        min_lp = np.amin(features, 0)
#        max_lp = np.amax(features, 0)
#        print(features)
#        print(min_lp, max_lp)
#        debug_terminate()
#        exit()
#        K = self.find_num_clusters()
#        centroids = self.select_initial_centroids()
#        old_centroids = np.zeros(centroids.shape)
#        candidates = []
#        errors = []
#        error = 0
#        m = 0
#        r = 0
#        print('round 0 centroids:', centroids)
#        # Kmeans main loop
#        while True:
#            m += 1
#            if m == 100 or m == -1 or (old_centroids == centroids).all():
#                if m == -1:
#                    print(red('empty cluster, skipping round'))
#                if m == 100:
#                    print(red('didn\'t converge, skipping round'))
#                    debug_breakpoint()
#                else:
#                    print(green('round', r), cyan('iteration', m))
#                    m = 0
#                    errors.append(error)
#                    candidates.append(centroids)
#                    centroids = self.select_initial_centroids(features, K)
#                    print('round', r, 'centroids:', centroids)
#                    old_centroids = np.zeros(centroids.shape)
#                r += 1
#                if r == 10:
#                    break
#            cluaters, error = self.cluster_data(features, centroids)
#            print('updated centroids:', centroids)
#            old_centroids = copy.deepcopy(centroids)
#            for k in range(K):
#                points = np.array([features[i] for i in range(len(features)) if clusters[i] == k])
#                print('cluster', k, 'points', points)
#                if len(points) == 0:
#                    m = -1
#                    print(red('invalid clustering'))
#                    continue
#                centroids[k] = np.mean(points, axis = 0)
#        choice = np.argmin(errors)
#        centroids = candidates[choice]
#        clusters, _ = self.cluster_data(features, centroids)
#        self.assign_genotypes(features, clusters, centroids)
#        self.plot_clusters(features, clusters)
#        return genotypes, centroids
#
#    def cluster_data(self, features, centroids):
#        error = 0
#        clusters = []
#        for i in range(len(features)):
#            distances = self.distance(features[i], centroids)
#            cluster = np.argmin(distances)
#            error += distances[cluster]
#            clusters.append(cluster)
#        return clusters, error
#
#    def assign_genotypes(self, features, clusters, centroids):
#        genotypes = []
#        norms = np.array([np.linalg.norm(centroid) for centroid in centroids])
#        print(norms)
#        centroids = []
#        for i in range(K):
#            centroid = np.argmin(norms)
#            norms[centroid] += 1
#            centroids.append(centroid)
#        print(centroids)
#        for i in range(len(features)):
#            if clusters[i] == centroids[0]:
#                genotypes.append('00')
#            elif clusters[i] == centroids[1]:
#                genotypes.append('10')
#            else:
#                genotypes.append('11')
#
#    def plot_clusters(self, features, clusters):
#        x = []
#        y = []
#        for i in range(len(features)):
#            x.append(clusters[i])
#            y.append(features[i][0])
#        try:
#            visualizer.violin(x, y, 'kmeans_' + track, self.get_current_job_directory(), 'cluster', 'value')
#        except Exception as e:
#            print(yellow(e))
#
#    def distance(self, a, b):
#        return np.linalg.norm(a - b, axis = 1)
#
#    def reduce(self):
#        files = {}
#        for path in self.paths:
#            name = path.split('/')[-1]
#            files[path] = open(os.path.join(self.get_current_job_directory(), name), 'w')
#        for batch in self.load_output():
#            for track in batch:
#                track = json.load(open(batch[track]))
#                for path in track:
#                    t = track[path]
#                    files[path].write(t['chrom'] + '\t' + str(t['begin']) + '\t' + str(t['end']) + '\t' + t['genotype'] + '\t' + t['id'] + '\t' + t['lp_value'] + '\n')
#
