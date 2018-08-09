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

from kmer import (
    bed,
    config,
    counter,
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import cplex
import numpy
import pybedtools

from Bio import pairwise2

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ContextSensitiveCountInnerKmersJob(counter.SimulationExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ContextSensitiveCountInnerKmersJob(job_name = 'CountInnerKmersJob_', previous_job_name = 'ExtractInnerKmersJob_', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        for i in range(0, self.num_threads):
            self.batch[i] = {}
        self.kmers = {}
        tracks = self.load_previous_job_results()
        index = {'inner_kmers': {}, 'novel_kmers': {}}
        self.sequence = {}
        for track in tracks:
            print(track)
            t = bed.track_from_name(track)
            chromosome = extract_chromosome(t.chrom)
            self.sequence[track] = chromosome[t.start - 1000 : t.end + 1000]
            print(len(self.sequence[track]))
            with open(tracks[track], 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers['inner_kmers']:
                    if not kmer in self.kmers:
                        self.kmers[kmer] = {track: 0}
                        self.kmers[reverse_complement(kmer)] = {track: 0}
                    else:
                        self.kmers[kmer][track] = 0
                        self.kmers[reverse_complement(kmer)][track] = 0
                for kmer in kmers['unique_inner_kmers']:
                    if not kmer in self.kmers:
                        self.kmers[kmer] = {track: 0}
                        self.kmers[reverse_complement(kmer)] = {track: 0}
                    else:
                        self.kmers[kmer][track] = 0
                        self.kmers[reverse_complement(kmer)][track] = 0
                with open(os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track + '.json'), 'w') as track_file:
                    json.dump(kmers, track_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            t = {}
            for track in tracks:
                t[track] = os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track + '.json')
            json.dump(t, json_file, indent = 4, sort_keys = True)
        chroms = range(1, 23)
        chroms.append('x')
        chroms.append('y')
        chroms = [2]
        chroms = ['chr' + str(chrom) for chrom in chroms]
        files = {}
        for chrom in chroms:
            for i in range (1, 3):
                for j in range(1, 3):
                    path = os.path.join(self.get_simulation_directory(), chrom + '_strand_' + str(i) + '.' + str(j) + '.fq')
                    files[path] = path
                    print(path)
        self.round_robin(files)

    def transform(self, track, track_name):
        c = config.Configuration()
        self.fastq_file = open(track, 'r')
        for read, name in self.parse_fastq():
            kmers = extract_kmers(c.ksize, read)
            scores = {}
            for kmer in kmers:
                if kmer in self.kmers: 
                    for track in self.kmers[kmer]:
                        if not track in scores:
                            print(track, len(self.sequence[track]))
                            alignments = pairwise2.align.globalxs(read, self.sequence[track], -1, -1)
                            gc.collect()
                            score = alignments[0][2]
                            score = 0
                            #print(score)
                            scores[track] = score
                        #for a in alignments:
                        #    print(format_alignment(*a))
                        if scores[track] > 0.9 * len(read):
                            self.kmers[kmer][track] += kmers[kmer]

    def merge_counts(self):
        c = config.Configuration()
        print('merging kmer counts ...')
        kmers = {}
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            if not os.path.isfile(path):
                print(red('couldn\'t find batch'), i, red('results will be unreliable'))
                continue
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    k = find_kmer(kmer, kmers)
                    if k:
                        for track in k:
                            kmers[k][track] += batch[kmer][track]
                    else:
                        kmers[kmer] = batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)
        return kmers

    def reduce(self):
        kmers = self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class ContextSensitiveIntegerProgrammingJob(map_reduce.BaseGenotypingJob):

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = ContextSensitiveIntegerProgrammingJob(job_name = 'IntegerProgramming_', previous_job_name = 'CountInnerKmersJob_' if c.simulation else 'ExtractInnerKmersJob_', category = 'programming', batch_file_prefix = 'unique_inner_kmers', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        with open(track, 'r') as json_file:
            kmers = json.load(json_file)
            inner_kmers = kmers['unique_inner_kmers']
            #inner_kmers.update(kmers['unique_inner_kmers'])
            if len(inner_kmers) == 0:
                print('no inner kmers found for', red(track_name))
                return None
            #inner_kmers = kmers['inner_kmers']
            for kmer in inner_kmers:
                if not kmer in self.inner_kmers:
                    self.inner_kmers[kmer] = {
                        'count': sum(map(lambda t: self.counts_provider.get_kmer_count(kmer)[t], self.counts_provider.get_kmer_count(kmer))),
                        'tracks': {},
                        'reference': self.reference_counts_provider.get_kmer_count(kmer)
                    }
                self.inner_kmers[kmer]['tracks'][track_name] = {'ref': inner_kmers[kmer]['track'], 'count': self.counts_provider.get_kmer_count(kmer)[track_name]}
            novel_kmers = {}
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner_kmers': {kmer: self.inner_kmers[kmer] for kmer in inner_kmers},
                    'novel_kmers': {kmer: self.novel_kmers[kmer] for kmer in novel_kmers},
                }, json_file, indent = 4, sort_keys = True)
        return path

    def incorporate_inner_kmers(self, problem):
        # the coverage of each event
        for index, track in enumerate(self.tracks):
            tokens = track['track'].split('_')
            problem.variables.add(names = ['c' + str(tokens[1])], ub = [1.0])
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.inner_kmers)],
            types = ['C'] * len(self.inner_kmers),
            ub = [kmer['count'] for kmer in self.inner_kmers],
            lb = [kmer['count'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks']) for kmer in self.inner_kmers])
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.inner_kmers)],
            obj = [1.0] * len(self.inner_kmers),
            types = ['C'] * len(self.inner_kmers))
        # constraints
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.inner_kmers):
            ind = list(map(lambda track: track, kmer['tracks'])) # C
            ind.append(len(self.tracks) + index) #E + str(index)
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
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = [len(self.tracks) + len(self.inner_kmers) + index, len(self.tracks) + index],
                    val = [1.0, 1.0],
                )],
                rhs = [0],
                senses = ['G']
            )
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = [len(self.tracks) + len(self.inner_kmers) + index, len(self.tracks) + index],
                    val = [1.0, -1.0],
                )],
                rhs = [0],
                senses = ['G']
            )
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(self.inner_kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
                print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
        return problem

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
