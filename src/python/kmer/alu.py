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

from kmer import (
    bed,
    config,
    gapped,
    counter,
    map_reduce,
    statistics,
    visualizer,
)

import pysam

from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractAluReadsJob(map_reduce.Job):

    _name = 'ExtractAluReadsJob'
    _category = 'programming'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractAluReadsJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        #extract_whole_genome()
        #self.chroms = extract_whole_genome()
        print(c.bam_file)
        self.bamfile = pysam.AlignmentFile(c.bam_file, "rb")
        self.tracks = bed.load_tracks_from_file_as_dict(c.bed_file)
        self.round_robin(self.tracks)
    
    def transform(self, track, track_name):
        gapped_kmers = {'outer': {}, 'inner': {}}
        slack = 70
        ref = ''
        reads = self.bamfile.fetch(region = track.chrom + ':' + str(track.begin - slack) + ':' + str(track.begin + slack))
        for read in reads:
            if read.query_alignment_length != len(read.query_sequence):
                seq = read.query_sequence
                if read.reference_start >= track.begin - slack and read.reference_start <= track.begin + slack:
                    cigar = read.cigartuples
                    clips = []
                    offset = 0
                    for c in cigar:
                        if c[0] == 4: #soft clip
                            clips.append((offset, offset + c[1]))
                        offset += c[1]
                    index = 0
                    for kmer in stream_kmers(37, False, seq):
                        if self.is_clipped((index, index + 37), clips):
                            if not kmer in gapped_kmers['outer']:
                                gapped_kmers['outer'][kmer] = 0
                            gapped_kmers['outer'][kmer] += 1
                        index += 1
        _gapped_kmers = {'outer': {}, 'inner': {}, 'ref': ref}
        for kmer in gapped_kmers['outer']:
            if gapped_kmers['outer'][kmer] >= 3:
                _gapped_kmers['outer'][kmer] = gapped_kmers['outer'][kmer]
        path = os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(_gapped_kmers, json_file, indent = 4)
        return path

    def is_clipped(self, kmer, clips):
        for clip in clips:
            if self.overlap(kmer, clip) >= 0 and self.overlap(kmer, clip) >= 10:
                return True
        return False

    def overlap(self, a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueAluGappedKmersJob(gapped.UniqueGappedKmersJob):

    _name = 'UniqueAluGappedKmersJob'
    _category = 'programming'
    _previous_job = ExtractAluReadsJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueAluGappedKmersJob(**kwargs)
        job.execute()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueAluGappedKmersScoringJob(gapped.UniqueGappedKmersScoringJob):

    _name = 'UniqueAluGappedKmersScoringJob'
    _category = 'programming'
    _previous_job = UniqueAluGappedKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueAluGappedKmersScoringJob(**kwargs)
        job.execute()

    def transform(self, sequence, chrom):
        c = config.Configuration()
        t = time.time()
        l = len(sequence)
        for index in range(0, l - 50):
            if index % 100000 == 1:
                s = time.time()
                p = index / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
            half_mer = sequence[index: index + c.hsize]
            if not half_mer in self.half_mers:
                continue
            for k in range(37, 38):
                kmer = canonicalize(sequence[index: index + k])
                kmer = kmer[:c.hsize] + kmer[-c.hsize:]
                if kmer in self.kmers:
                    self.kmers[kmer]['count'] += 1

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class SelectUniqueAluGappedKmersJob(gapped.SelectUniqueGappedKmersJob):

    _name = 'SelectUniqueAluGappedKmersJob'
    _category = 'programming'
    _previous_job = UniqueAluGappedKmersScoringJob
    _counter_mode = 2

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = SelectUniqueAluGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.half_mers = {}
        tracks = self.load_previous_job_results()
        m = {}
        for track in tracks:
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['outer']:
                    if kmers['outer'][kmer]['count'] == 0:
                        left = kmer[:c.hsize]
                        right = kmer[-c.hsize:]
                        self.kmers[kmer] = {'tracks': kmers['outer'][kmer]['tracks'], 'side': 'outer', 'count': {}}
                        for i in range(0, 11):
                            self.kmers[kmer]['count'][i] = 0
                        if not left in self.half_mers:
                            self.half_mers[left] = {}
                        self.half_mers[left][right] = kmer 
                        left = reverse_complement(left)
                        right = reverse_complement(right)
                        if not right in self.half_mers:
                            self.half_mers[right] = {}
                        self.half_mers[right][left] = kmer 
                        if track not in m:
                            m[track] = 0
                        m[track] += 1
        visualizer.histogram(list(map(lambda track: m[track], m)), 'gapped_kmers_per_alu', self.get_current_job_directory(), 'number of kmers', 'number of events')
        print(len(tracks), 'total events')
        print(len(self.kmers), 'kmers')
        print(len(m), 'events with kmers')
        self.export_accelerator_input()
        self.round_robin()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountUniqueAluGappedKmersJob(gapped.CountUniqueGappedKmersJob):

    _name = 'CountUniqueAluGappedKmersJob'
    _category = 'programming'
    _previous_job = SelectUniqueAluGappedKmersJob
    _counter_mode = 1

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountUniqueAluGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        tracks = self.load_previous_job_results()
        self.half_mers = {}
        n = 0
        for track in tracks:
            n += 1
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    for kmer in kmers[side]:
                        if kmers[side][kmer]['gap'] == 5:
                            left = kmer[:c.hsize]
                            right = kmer[-c.hsize:]
                            self.kmers[kmer] = kmers[side][kmer]
                            self.kmers[kmer]['count'] = 0
                            self.kmers[kmer]['doubt'] = 0
                            if not left in self.half_mers:
                                self.half_mers[left] = {}
                            self.half_mers[left][right] = kmer 
                            left = reverse_complement(left)
                            right = reverse_complement(right)
                            if not right in self.half_mers:
                                self.half_mers[right] = {}
                            self.half_mers[right][left] = kmer
        print(len(self.kmers), 'kmers')
        self.export_accelerator_input()
        self.round_robin()

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class AluGappedKmersIntegerProgrammingJob(gapped.GappedKmersIntegerProgrammingJob):

    _name = 'AluGappedKmersIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = CountUniqueAluGappedKmersJob
    _kmer_type = 'gapped'

    @staticmethod
    def launch(**kwargs):
        job = AluGappedKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def generate_linear_program(self):
        print('generating linear program')
        c = config.Configuration()
        globals()['cplex'] = __import__('cplex')
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
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            lb = [(kmer['count'] - kmer['coverage']) for kmer in self.lp_kmers],
        )
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [1.0] * len(self.lp_kmers),
        )
        # constraints
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.lp_kmers):
            ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
            ind.append(len(self.tracks) + index)
            val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
            val.append(1.0)
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [kmer['count']],
                senses = ['E']
            )
            self.add_error_absolute_value_constraints(problem, index)
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(self.lp_kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
        return problem

    def round_genotype(self, c):
        if c > 0.75:
            return (1.0, '11')
        elif c > 0.25:
            return (0.5, '10')
        else:
            return (0.0, '00')

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class AluGappedKmersFractionalProgrammingJob(gapped.GappedKmersIntegerProgrammingJob):

    _name = 'AluGappedKmersFractionalProgrammingJob'
    _category = 'programming'
    _previous_job = CountUniqueAluGappedKmersJob
    _kmer_type = 'gapped'

    @staticmethod
    def launch(**kwargs):
        job = AluGappedKmersFractionalProgrammingJob(**kwargs)
        job.execute()

    def solve(self):
        helper = AluGappedKmersFractionalProgrammingJob.FractionalProgramHelper()
        helper.lp_kmers = self.lp_kmers
        helper.tracks = self.tracks
        helper.execute()
        self.export_solution()

    def export_solution(self):
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            for track in self.tracks:
                with open(os.path.join(self.get_current_job_directory(), 'solution_' + track + '.json'), 'r') as json_file:
                    solution = json.load(json_file)['variables']
                t = bed.track_from_name(track)
                c = solution[0]
                w = solution[-1]
                print(c)
                #g = self.round_genotype(c / w)
                g = self.round_genotype(c / 2.0)
                w = 1.0
                bed_file.write(t.chrom + '\t' +
                    str(t.begin) + '\t' +
                    str(t.end) + '\t' +
                    str(g[1]) + '\t' +
                    str(solution[1]) + '\t' +
                    str(solution[2]) + '\t' +
                    str(solution[3]) + '\n')
                    #str(1.0 / w) + '\t' +
                    #str(c / w) + '\t' +
                    #str(len(solution) / 2) + '\n')
                with open(os.path.join(self.get_current_job_directory(), 'solution_' + track + '.json'), 'w') as json_file:
                    json.dump({'variables': [s / w if s != w else 1.0 / s for s in solution]}, json_file, indent = 4, sort_keys = True)

    def round_genotype(self, c):
        if c > 0.75:
            return (1.0, '11')
        elif c > 0.25:
            return (0.5, '10')
        else:
            return (0.0, '00')

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    class FractionalProgramHelper(gapped.GappedKmersIntegerProgrammingJob):

        _name = 'AluGappedKmersFractionalProgrammingJob'
        _category = 'programming'
        _previous_job = CountUniqueAluGappedKmersJob
        _kmer_type = 'gapped'

        def load_inputs(self):
            self.round_robin(self.tracks)

        def transform(self, _, track):
            print('solving for', cyan(track))
            #problem = self.generate_linear_program(track)
            #print('program_' + track + '.lp')
            #problem.write(os.path.join(self.get_current_job_directory(), 'program_' + str(track) + '.lp'))
            #problem.solve()
            #solution = problem.solution.get_values()
            solution = self.estimate_maximum_likelihood(track)
            with open(os.path.join(self.get_current_job_directory(), 'solution_' + track + '.json'), 'w') as json_file:
                json.dump({'variables': solution}, json_file, indent = 4, sort_keys = True)
            return None

        def reduce(self):
            pass

        def output_batch(self, batch):
            exit()

        def estimate_maximum_likelihood(self, track):
            c = config.Configuration()
            lp_kmers = list(filter(lambda kmer: track in kmer['tracks'], self.lp_kmers))
            distributions = [statistics.NormalDistribution(0, 3), statistics.NormalDistribution(c.coverage / 2, 8), statistics.NormalDistribution(c.coverage, 15)]
            likelihoods = [0, 0, 0]
            for index, kmer in enumerate(lp_kmers):
                for i, distribution in enumerate(distributions):
                    likelihoods[i] += distribution.log_pmf(kmer['count'])
            index, value = max(enumerate(likelihoods), key = operator.itemgetter(1))
            return [index] + likelihoods

        def generate_linear_program(self, track):
            print('generating linear program')
            c = config.Configuration()
            globals()['cplex'] = __import__('cplex')
            problem = cplex.Cplex()
            problem.objective.set_sense(problem.objective.sense.maximize)
            lp_kmers = list(filter(lambda kmer: track in kmer['tracks'], self.lp_kmers))
            # the coverage of each event
            tokens = track.split('_')
            problem.variables.add(names = ['y_c' + tokens[1]],
                lb = [-100000],
            )
            # the real-valued error parameter for inner_kmer
            problem.variables.add(names = ['y_e' + str(index) for index, kmer in enumerate(lp_kmers)],
                lb = [-100000 for index, kmer in enumerate(lp_kmers)],
            )
            # absolute value of the inner_kmer error parameter
            problem.variables.add(names = ['y_l' + str(index) for index, kmer in enumerate(lp_kmers)],
                obj = [-1.0] * len(lp_kmers),
            )
            # t variable
            problem.variables.add(names = ['t'], lb = [-100000])
            # constraints
            n = 0
            start = time.time()
            ind = [0]
            ind.append(1 + 2 * len(lp_kmers))
            #ind.append(1 + index)
            val = [17, 0.01]
            #val.append([kmer['count']] * (1 + 2 * len(lp_kmers)))
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [1],
                senses = ['E']
            )
            for index, kmer in enumerate(lp_kmers):
                ind = [0]
                ind.append(1 + index)
                ind.append(1 + 2 * len(lp_kmers))
                val = [kmer['coverage']]
                val.append(1.0)
                val.append(-1 * kmer['count'])
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [0],
                    senses = ['E']
                )
                # error absolute value
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = [1 + len(lp_kmers) + index, 1 + index],
                        val = [-1.0, 1.0],
                    )],
                    rhs = [0],
                    senses = ['L']
                )
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = [1 + len(lp_kmers) + index, 1 + index],
                        val = [-1.0, -1.0],
                    )],
                    rhs = [0],
                    senses = ['L']
                )
                n = n + 1
                if n % 1000 == 0:
                    t = time.time()
                    p = float(n) / len(self.lp_kmers)
                    eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
            return problem

