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
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from nebula.kmers import *
from nebula.logger import *
from nebula.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class VertexCover1000Genomes(map_reduce.Job):

    _name = 'VertexCover1000Genomes'
    _category = 'misc'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = VertexCover1000Genomes(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        vcf_path = c.input
        self.tracks = []
        self.genomes = {}
        t_index = 0
        with open(vcf_path, 'r') as vcf_file:
            header = vcf_file.readline()
            tokens = header.split()
            for index, genome in enumerate(tokens[9:]):
                self.genomes[index] = {'name': genome, 'tracks': {}}
            tracks = vcf_file.readlines()
            n = 0
            m = 0
            for track in tracks:
                n += 1
                print(n, 'out of', len(tracks))
                try:
                    tokens = track.split()
                    fields = {p[0]: p[1] for p in list(map(lambda x: (x[:x.find('=')], x[x.find('=') + 1:]), tokens[7].split(';')))}
                    if fields['SVTYPE'] != 'DEL':
                        continue
                    t = bed.BedTrack(chrom = 'chr' + tokens[0], begin = tokens[1], end = fields['END'], genomes = list(filter(lambda i: tokens[i + 9] != '0|0', range(0, len(tokens) - 9))), index = t_index)
                    t_index += 1
                    if len(t.genomes) > 50:
                        for genome in t.genomes:
                            self.genomes[genome]['tracks'][t.index] = True
                        self.tracks.append(t)
                        m += 1
                        print(m)
                except Exception as e:
                    print(red(e))
                    print(track)
                    debug_breakpoint()
        print('Tracks with AC > 5%:', len(self.tracks))
        #self.generate_lp()
        #debug_breakpoint()
        self.solve_greedy()
        exit()

    def solve_greedy(self):
        print('solving greedy')
        targets = []
        covered = {}
        print('Trying to cover', len(self.tracks), 'tracks')
        while len(covered) != len(self.tracks):
            m = -1
            select = None
            s = 0
            for genome in self.genomes:
                s = len(list(filter(lambda x: x not in covered, self.genomes[genome]['tracks'])))
                if s > m:
                    m = s
                    select = genome
            print('added', self.genomes[select]['name'], 'covering', m, 'new tracks')
            for track in self.genomes[select]['tracks']:
                covered[track] = True
            targets.append(select)
            self.genomes.pop(select, None)
            print('Coverage:', len(covered))
        print('Minimum number of genomes:', len(targets), 'out of', len(self.genomes))

    def generate_lp(self):
        globals()['cplex'] = __import__('cplex')
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        for genome in self.genomes:
            problem.variables.add(names = [self.genomes[genome]['name']], ub = [1], lb = [0], types = [problem.variables.type.integer], obj = [1.0])
        targets = []
        for index, track in enumerate(self.tracks):
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                        ind = track.genomes,#list(map(lambda g: self.genomes[g]['index'], track.genomes)),
                        val = [1] * len(track.genomes)
                    )],
                    rhs = [1],
                    senses = ['G']
                )
        problem.write(os.path.join(self.get_current_job_directory(), 'program.lp'))
        problem.solve()
        solutions = problem.solution.get_values()
        targets = []
        for index, sol in enumerate(solutions):
            if int(sol) == 1:
                targets.append(self.genomes[index])
        #json_print(targets)
        print('Minimum number of genomes:', len(targets), 'out of', len(self.genomes))
        #with open(os.path.join(self.get_current_job_directory(), 'solution.json'), 'w') as json_file:
        #    json.dump(targets, json_file, indent = 4)







