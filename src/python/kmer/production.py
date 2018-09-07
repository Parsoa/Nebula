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

class MixIntegerProgrammingJob(programming.IntegerProgrammingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = MixIntegerProgrammingJob(job_name = 'MixGappedInnerKmersIntegerProgrammingJob_', previous_job_name = '', category = 'programming', batch_file_prefix = 'batch',
            include_gapped = True,
            include_unique_inner = True,
            include_non_unique_inner = True,
        )
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare(self):
        self.job_name = 'Mix'
        if self.include_gapped:
            self.job_name += 'Gapped'
        if self.include_unique_inner:
            self.job_name += 'UniqueInner'
        if self.include_non_unique_inner:
            self.job_name += 'NoneUniqueInner'
        self.job_name += 'KmersIntegerProgrammingJob_'
        print(self.job_name)

    def execute(self):
        c = config.Configuration()
        self.prepare()
        self.create_output_directories()
        self.bedtools = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'all.bed'))], key = lambda track: track.start)
        self.all_tracks = {}
        n = 0
        for b in self.bedtools:
            track = re.sub(r'\s+', '_', str(b).strip()).strip()
            print(track)
            self.all_tracks[track] = n
            n += 1
        self.gapped_kmers = []
        self.gapped_tracks = {}
        self.unique_inner_kmers = []
        self.unique_inner_tracks = {}
        self.non_unique_inner_kmers = []
        self.non_unique_inner_tracks = {}
        #
        self.load_gapped_kmers()
        self.load_unique_inner_kmers()
        self.load_non_unique_inner_kmers()
        self.merge_kmers()
        self.solve()
        exit()

    def load_gapped_kmers(self):
        if not self.include_gapped:
            return
        print(cyan('=============================================================================================='))
        print(cyan('Loading gapped kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'r') as json_file:
                payload = json.load(json_file)
                self.gapped_kmers = payload['gapped_kmers']
                self.gapped_tracks = payload['tracks']
            return
        self.gapped_kmers_solver = gapped.GappedKmersIntegerProgrammingJob(job_name = 'GappedKmersIntegerProgrammingJob_', previous_job_name = 'CountUniqueGappedKmersJob_', category = 'programming', batch_file_prefix = 'gapped_kmers', kmer_type = 'gapped')
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

    def load_unique_inner_kmers(self):
        if not self.include_unique_inner:
            return
        print(cyan('=============================================================================================='))
        print(cyan('Loading unique inner kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'unique_inner_kmers.json'), 'r') as json_file:
                payload = json.load(json_file)
                self.unique_inner_kmers = payload['unique_inner_kmers']
                self.unique_inner_tracks = payload['tracks']
            return
        self.unique_inner_kmers_solver = programming.IntegerProgrammingJob(job_name = 'IntegerProgrammingJob_', previous_job_name = 'CountInnerKmersJob_' if c.simulation else 'ExtractInnerKmersJob_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', kmer_type = 'unique_inner')
        self.unique_inner_kmers_solver.create_output_directories()
        self.unique_inner_kmers_solver.load_inputs()
        self.unique_inner_kmers_solver.distribute_workload()
        self.unique_inner_kmers_solver.wait_for_children()
        self.unique_inner_kmers = copy.deepcopy(self.unique_inner_kmers_solver.index_kmers())
        self.unique_inner_tracks = copy.deepcopy(self.unique_inner_kmers_solver.index_tracks())
        with open(os.path.join(self.get_current_job_directory(), 'unique_inner_kmers.json'), 'w') as json_file:
            json.dump({'unique_inner_kmers': self.unique_inner_kmers, 'tracks': self.unique_inner_tracks}, json_file, indent = 4)
        self.unique_inner_kmers_solver.calculate_residual_coverage()
        self.unique_inner_kmers_solver.solve()
 
    def load_non_unique_inner_kmers(self):
        if not self.include_non_unique_inner:
            return False
        print(cyan('=============================================================================================='))
        print(cyan('Loading inner kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'non_unique_inner_kmers.json'), 'r') as json_file:
                payload = json.load(json_file)
                self.non_unique_inner_kmers = payload['non_unique_inner_kmers']
                self.non_unique_inner_tracks = payload['tracks']
            return
        self.non_unique_inner_kmers_solver = reduction.LociIndicatorKmersIntegerProgrammingJob(job_name = 'LociIndicatorKmersIntegerProgrammingJob_', previous_job_name = 'CountLociIndicatorKmersJob_', category = 'programming', batch_file_prefix = 'indicator_kmers', kmer_type = 'non_unique_inner')
        self.non_unique_inner_kmers_solver.create_output_directories()
        self.non_unique_inner_kmers_solver.load_inputs()
        self.non_unique_inner_kmers_solver.distribute_workload()
        self.non_unique_inner_kmers_solver.wait_for_children()
        self.non_unique_inner_kmers = copy.deepcopy(self.non_unique_inner_kmers_solver.index_kmers())
        self.non_unique_inner_tracks = copy.deepcopy(self.non_unique_inner_kmers_solver.index_tracks())
        with open(os.path.join(self.get_current_job_directory(), 'non_unique_inner_kmers.json'), 'w') as json_file:
            json.dump({'non_unique_inner_kmers': self.non_unique_inner_kmers, 'tracks': self.non_unique_inner_tracks}, json_file, indent = 4)
        self.non_unique_inner_kmers_solver.calculate_residual_coverage()
        self.non_unique_inner_kmers_solver.solve()

    def merge_kmers(self):
        self.tracks = {}
        self.tracks.update(self.non_unique_inner_tracks)
        self.tracks.update(self.unique_inner_tracks)
        self.tracks.update(self.gapped_tracks)
        #self.find_tracks_with_no_signal()
        #for track in self.tracks.keys():
        #    if track in self.unique_inner_tracks or track in self.gapped_tracks:
        #        self.tracks.pop(track, None)
        n = 0
        tmp = sorted([t for t in self.tracks])
        for track in tmp:
            self.tracks[track] = {'index': n, 'unique_inner': 0, 'gapped': 0, 'non_unique_inner': 0}
            n += 1
        self.lp_kmers = self.unique_inner_kmers + self.gapped_kmers + self.non_unique_inner_kmers
        #self.lp_kmers = self.non_unique_inner_kmers
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
            kmer['coverage'] = 42 if kmer['type'] == 'gapped' else 50

    def generate_linear_program(self):
        c = config.Configuration()
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        for track in self.tracks:
            tokens = track.split('_')
            problem.variables.add(names = ['c' + str(tokens[1])],
                ub = [1.0],
            )
        # the real-valued error parameter for kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            lb = [(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in self.lp_kmers]
        )
        # absolute value of the kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [1.0 / kmer['reference'] for index, kmer in enumerate(self.lp_kmers)]
        )
        m = 0
        n = 0
        for index, kmer in enumerate(self.lp_kmers):
            self.add_error_absolute_value_constraints(problem, index)
            #if kmer['type'] == 'non_unique_inner':
            #    if len(filter(lambda track: track not in self.tracks, kmer['tracks'])) != 0:
            #        problem.linear_constraints.add(
            #            lin_expr = [cplex.SparsePair(
            #                ind = [len(self.tracks) + index],
            #                val = [1.0],
            #            )],
            #            rhs = [0],
            #            senses = ['E']
            #        )
            #        problem.variables.set_lower_bounds(len(self.tracks) + index, 0)
            #        continue
            if kmer['type'] == 'non_unique_inner':
                if len(kmer['tracks']) == 1:
                    track = kmer['tracks'].keys()[0]
                    if self.tracks[track]['gapped'] + self.tracks[track]['unique_inner'] > 5:
                        m += 1
                        problem.linear_constraints.add(
                            lin_expr = [cplex.SparsePair(
                                ind = [len(self.tracks) + index],
                                val = [1.0],
                            )],
                            rhs = [0],
                            senses = ['E']
                        )
                        problem.variables.set_lower_bounds(len(self.tracks) + index, 0)
                        continue
                    else:
                        n += 1
                        print(yellow(track))
            for track in kmer['tracks']:
                self.tracks[track][kmer['type']] += 1
            if kmer['type'] == 'gapped' and kmer['side'] == 'outer':
                # (1 - T)xR + E = C -> -TxR + E = C - R
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - kmer['coverage'] * kmer['residue'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['E']
                )
            else:
                # TxR + E = C
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
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
        print(m, n)
        return problem

    def solve(self):
        programming.IntegerProgrammingJob.solve(self)
        for track in self.tracks:
            self.tracks[track] = {'unique_inner_kmers': {}, 'non_unique_inner_kmers': {}, 'gapped_kmers': {}}
        for kmer in self.unique_inner_kmers:
            for track in kmer['tracks']:
                self.tracks[track]['unique_inner_kmers'][kmer['kmer']] = kmer
        for kmer in self.non_unique_inner_kmers:
            for track in kmer['tracks']:
                self.tracks[track]['non_unique_inner_kmers'][kmer['kmer']] = kmer
        for kmer in self.gapped_kmers:
            for track in kmer['tracks']:
                self.tracks[track]['gapped_kmers'][kmer['kmer']] = kmer
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
