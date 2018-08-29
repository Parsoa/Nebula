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

class MixGappedInnerKmersIntegerProgrammingJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = MixGappedInnerKmersIntegerProgrammingJob(job_name = 'MixGappedInnerKmersIntegerProgrammingJob_', previous_job_name = '', category = 'programming', batch_file_prefix = 'batch', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def execute(self):
        c = config.Configuration()
        self.bedtools = sorted([track for track in pybedtools.BedTool(c.bed_file)], key = lambda track: track.start)
        self.tracks = {}
        n = 0
        for b in self.bedtools:
            track = re.sub(r'\s+', '_', str(b).strip()).strip()
            print(track)
            self.tracks[track] = n
            n += 1
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
                self.gapped_kmers = json.load(json_file)['gapped_kmers']
                return
        self.gapped_kmers_solver = gapped.GappedKmersIntegerProgrammingJob(job_name = 'GappedKmersIntegerProgrammingJob_', previous_job_name = 'CountUniqueGappedKmersJob_', category = 'programming', batch_file_prefix = 'unique_gapped_kmers')
        self.gapped_kmers_solver.create_output_directories()
        self.gapped_kmers_solver.load_inputs()
        self.gapped_kmers_solver.distribute_workload()
        self.gapped_kmers_solver.wait_for_children()
        self.gapped_kmers = copy.deepcopy(self.gapped_kmers_solver.index_kmers())
        with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
            json.dump({'gapped_kmers': self.gapped_kmers}, json_file, indent = 4)
        self.gapped_kmers_solver.calculate_residual_coverage()
        self.gapped_kmers_solver.solve()

    def load_inner_kmers(self):
        print(cyan('=============================================================================================='))
        print(cyan('Loading inner kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'r') as json_file:
                self.inner_kmers = json.load(json_file)['inner_kmers']
            return
        self.inner_kmers_solver = programming.IntegerProgrammingJob(job_name = 'IntegerProgrammingJob_', previous_job_name = 'CountInnerKmersJob_' if c.simulation else 'ExtractInnerKmersJob_', category = 'programming', batch_file_prefix = 'unique_inner_kmers')
        self.inner_kmers_solver.create_output_directories()
        self.inner_kmers_solver.load_inputs()
        self.inner_kmers_solver.distribute_workload()
        self.inner_kmers_solver.wait_for_children()
        self.inner_kmers = copy.deepcopy(self.inner_kmers_solver.index_kmers())
        with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
            json.dump({'inner_kmers': self.inner_kmers}, json_file, indent = 4)
        self.inner_kmers_solver.calculate_residual_coverage()
        self.inner_kmers_solver.solve()

    def merge_kmers(self):
        self.kmers = self.inner_kmers + self.gapped_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump({'kmers': self.kmers}, json_file, indent = 4)
        #print(self.kmers)
        print(len(self.kmers))
        print(cyan('=============================================================================================='))
        print(cyan('Solving composite model on'), blue(len(self.kmers)), cyan('kmers and '), blue(len(self.tracks)), cyan('tracks'))
        print(cyan('=============================================================================================='))
        for kmer in self.kmers:
            r = 0
            print(kmer['tracks'])
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['residue'] = 0 if kmer['type'] == 'gapped' else kmer['reference'] - r
            kmer['coverage'] = c.coverage

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
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.kmers)],
            lb = [(kmer['count'] - c.coverage * kmer['residue'] - c.coverage * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in self.kmers]
        )
        # absolute value of the kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.kmers)],
            obj = [1.0] * len(self.kmers),
        )
        for index, kmer in enumerate(self.kmers):
            if kmer['type'] == 'gapped' and kmer['side'] == 'outer':
                # (1 - T)xR + E = C -> -TxR + E = C - R
                ind = list(map(lambda track: self.tracks[track], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: -1 * c.coverage * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - c.coverage * kmer['residue'] - sum(list(map(lambda track: c.coverage * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['E']
                )
            else:
                # TxR + E = C
                ind = list(map(lambda track: self.tracks[track], kmer['tracks']))
                ind.append(len(self.tracks) + index)
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
        return problem

    def solve(self):
        problem = self.generate_linear_program()
        problem.write(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_program.lp'))
        problem.solve()
        solution = problem.solution.get_values()
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_solution.json'), 'w') as json_file:
            json.dump({'variables': problem.solution.get_values()}, json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_merge.bed'), 'w') as bed_file:
            for track in self.tracks:
                tokens = track.split('_')
                s = int(round(2 * solution[self.tracks[track]]))
                s = '(0, 0)' if s == 2 else '(1, 0)' if s == 1 else '(1, 1)'
                bed_file.write(tokens[0] + '\t' + #0
                            tokens[1] + '\t' + #1
                            tokens[2] + '\t' + #2
                            s + '\t' + #3
                            str(solution[self.tracks[track]]) + '\t' + #4
                            #str(len(track['inner_kmers'])) + '\t' + #5
                            self.batch_file_prefix + '\n') #6

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
