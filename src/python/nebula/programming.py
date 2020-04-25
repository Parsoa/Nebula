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
from nebula.logger import *
from nebula.chromosomes import *
print = pretty_print

from pulp import *

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class IntegerProgrammingJob(map_reduce.Job):

    _name = 'IntegerProgrammingJob'
    _category = 'preprocessing'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = IntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.round_robin(self.tracks)
        self.resume_from_reduce = False
        self.lp_kmers = {}

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.lp_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def reduce(self):
        c = config.Configuration()
        self.index_kmers()
        self.index_tracks()
        self.calculate_residual_coverage()
        self.solve()

    def index_kmers(self):
        c = config.Configuration()
        self.tracks = {}
        self.lp_kmers = []
        index = {}
        for batch in self.load_output():
            kmers = batch
            for kmer in kmers:
                if not kmer in index:
                    index[kmer] = len(self.lp_kmers)
                    self.lp_kmers.append(copy.deepcopy(kmers[kmer]))
                    self.lp_kmers[len(self.lp_kmers) - 1]['kmer'] = kmer
                    for track in kmers[kmer]['tracks']:
                        if not track in self.tracks:
                            self.tracks[track] = bed.track_from_id(track)
                            if kmers[kmer]['type'] == 'inner':
                                self.tracks[track]['confidence'] = 'LOW'
                            else:
                                self.tracks[track]['confidence'] = 'HIGH'
        return self.lp_kmers

    def index_tracks(self):
        n = 0
        tmp = sorted([t for t in self.tracks])
        for track in tmp:
            self.tracks[track]['index'] = n
            self.tracks[track]['kmers'] = []
            n += 1
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['kmers'].append(index)
        return self.tracks

    def add_mps_error_absolute_value_constraints(self, problem, variables, index):
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], 1.0), (variables[len(self.tracks) + index], 1.0)])
        problem += LpConstraint(expr, LpConstraintGE, 'abs_1_' + str(index), 0) 
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], 1.0), (variables[len(self.tracks) + index], -1.0)])
        problem += LpConstraint(expr, LpConstraintGE, 'abs_2_' + str(index), 0) 

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
        self.round_lp()
        self.export_solution()
        self.export_kmers()
        if c.rum:
            self.verify_genotypes()
        return self.tracks

    # COIN doesn't supply values for certain variables
    def import_lp_values(self, path = 'solution.mps'):
        c = config.Configuration()
        self.solution = [0.0] * (len(self.tracks) + 2 * len(self.lp_kmers))
        var_index = {}
        regex = re.compile('[^a-zA-Z0-9]')
        for track in self.tracks:
            name = regex.sub('_', track)
            var_index[name] = self.tracks[track]['index']
        with open(os.path.join(self.get_current_job_directory(), path), 'r') as f:
            status = f.readline()
            objective = f.readline()
            line = f.readline()
            while(line):
                tokens = line.split()
                name = tokens[1]
                index = int(tokens[0])
                value = float(tokens[2])
                # tracks are indexed lexicographically, different from iteration index
                if name[0] == 'c':
                    self.solution[var_index[name[1:]]] = value
                else:
                    self.solution[index] = value
                line = f.readline()

    def round_lp(self):
        c = config.Configuration()
        print('Rounding', len(self.tracks), 'tracks')
        for track in self.tracks:
            index = self.tracks[track]['index']
            self.tracks[track]['lp_kmers'] = []
            self.tracks[track]['lp_value'] = self.solution[index]
            self.tracks[track]['lp_genotype'] = self.round_genotype(self.solution[index], self.tracks[track]['svtype'])[1]

    def export_solution(self):
        c = config.Configuration()
        name = 'merge.bed' if not c.cgc else 'genotypes_' + (c.fastq.split('/')[-1] if c.fastq else c.bam.split('/')[-1]) + '.bed'
        path = os.path.join(c.workdir if c.cgc else self.get_current_job_directory(), name)
        with open(path, 'w') as bed_file:
            bed_file.write('#CHROM\tBEGIN\tEND\tID\tSVTYPE\tLP\tCONFIDENCE\tGENOTYPE\n')
            for t in bed.sort_tracks(self.tracks):
                bed_file.write('\t'.join([str(x) for x in [t.chrom, t.begin, t.end, t.id, t.svtype, t['lp_value'], t['confidence'], t['lp_genotype']]]) + '\n')

    def export_kmers(self):
        c = config.Configuration()
        tracks = {}
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                if not track in tracks:
                    tracks[track] = {}
                tracks[track][kmer['kmer']] = kmer
                self.tracks[track]['lp_kmers'].append(kmer)
        for track in tracks:
            path = os.path.join(self.get_current_job_directory(), track + '.json')
            with open(path, 'w') as json_file:
                json.dump(tracks[track], json_file, indent = 4)
        if c.cgc:
            name = 'lp_kmer_' + (c.fastq.split('/')[-1] if c.fastq else c.bam.split('/')[-1]) + '.json'
            with open(os.path.join(c.workdir, name), 'w') as json_file:
                json.dump({kmer['kmer']: kmer for kmer in self.lp_kmers}, json_file, indent = 4, sort_keys = True)
