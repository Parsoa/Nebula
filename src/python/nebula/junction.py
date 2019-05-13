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

from itertools import chain

from nebula import (
    bed,
    config,
    gapped,
    counter,
    reduction,
    map_reduce,
    statistics,
    visualizer,
)

import pysam

from nebula.kmers import *
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractJunctionKmersJob(map_reduce.Job):

    _name = 'ExtractJunctionKmersJob'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractJunctionKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def execute(self):
        c = config.Configuration()
        self.pre_process()
        self.create_output_directories()
        self.find_thread_count()
        self.prepare()
        self.load_inputs()

    def load_inputs(self):
        c = config.Configuration()
        self.contigs = pysam.AlignmentFile(c.contigs, "rb")
        self.alignments = pysam.AlignmentFile(c.bam, "rb")
        self.tracks = c.tracks
        self.extract_kmers()
        self.num_threads = 1

    def extract_kmers(self):
        c = config.Configuration()
        contig_index = {}
        output = {}
        for track in self.tracks:
            contig_index[self.tracks[track].contig] = self.tracks[track]
        n = 0
        print(len(contig_index), 'total tracks')
        for read in self.contigs.fetch():
            if read.query_name in contig_index:
                t = contig_index[read.query_name]
                kmers = {}
                kmers.update(self.extract_assembly_kmers(t, read.query_sequence))
                kmers.update(self.extract_mapping_kmers(t))
                if len(kmers) > 0:
                    path = os.path.join(self.get_current_job_directory(), t.id + '.json')
                    with open(path, 'w') as json_file:
                        json.dump(kmers, json_file, indent = 4)
                    output[t.id] = path
                    n += 1
                if n % 1000 == 0:
                    print(n)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump(output, json_file, indent = 4)

    def extract_assembly_kmers(self, track, assembly):
        c = config.Configuration()
        junction_kmers = {}
        if track.svtype == 'DEL':
            l = list(range(int(track.contig_start) - 24, int(track.contig_start) - 8))
        if track.svtype == 'INS':
            l = list(range(int(track.contig_start) - 24, int(track.contig_start) - 8)) + list(range(int(track.contig_end) - 24, int(track.contig_end) - 8))
        for i in l:
            kmer = canonicalize(assembly[i: i + c.ksize])
            if not kmer in junction_kmers:
                junction_kmers[kmer] = {
                    'count': 0,
                    'source': 'assembly',
                    'loci': {
                        'junction_' + track.id: {
                            'masks': {
                                'left': assembly[i - c.ksize: i],
                                'right': assembly[i + c.ksize: i + c.ksize + c.ksize]
                            }
                        }
                    }
                }
            junction_kmers[kmer]['count'] += 1
        junction_kmers = {kmer: junction_kmers[kmer] for kmer in junction_kmers if junction_kmers[kmer]['count'] == 1}
        return junction_kmers

    def extract_mapping_kmers(self, track):
        c = config.Configuration()
        junction_kmers = {}
        slack = 100
        stride = c.ksize
        #if not c.simulation:
        #    track = track.lift()
        if track.svtype == 'DEL':
            reads = chain(self.alignments.fetch(track.chrom, track.begin - slack, track.begin + slack), self.alignments.fetch(track.chrom, track.end - slack, track.end + slack))
        if track.svtype == 'INS':
            reads = self.alignments.fetch(track.chrom, track.begin - slack, track.begin + slack)
        n = 0
        strid = c.ksize
        for read in reads:
            n += 1
            if read.query_alignment_length != len(read.query_sequence): # not everything was mapped
                seq = read.query_sequence
                if read.reference_start >= track.begin - slack and read.reference_start <= track.end + slack:
                    cigar = read.cigartuples
                    clips = []
                    offset = 0
                    for c in cigar:
                        if c[0] == 4: #soft clip
                            clips.append((offset, offset + c[1]))
                        offset += c[1]
                    index = 0
                    for kmer in stream_kmers(stride, False, True, seq):
                        if self.is_clipped((index, index + stride), clips):
                            if not kmer in junction_kmers:
                                junction_kmers[kmer] = {
                                    'count': 0,
                                    'source': 'mapping',
                                    'loci': {
                                        'junction_' + track.id: {
                                            'masks': {
                                                'left': None,
                                                'right': None
                                            }
                                        }
                                    }
                                }
                            if len(seq[index - stride: index]) == stride:
                                junction_kmers[kmer]['loci']['junction_' + track.id]['masks']['left'] = seq[index - stride: index]
                            if len(seq[index + stride: index + stride + stride]) == stride:
                                junction_kmers[kmer]['loci']['junction_' + track.id]['masks']['right'] = seq[index + stride: index + stride + stride]
                            junction_kmers[kmer]['count'] += 1
                        index += 1
        _junction_kmers = {} 
        for kmer in junction_kmers:
            if junction_kmers[kmer]['count'] >= 3 and junction_kmers[kmer]['loci']['junction_' + track.id]['masks']['right'] and junction_kmers[kmer]['loci']['junction_' + track.id]['masks']['left']:
                _junction_kmers[kmer] = junction_kmers[kmer]
        return _junction_kmers

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

class UniqueJunctionKmersJob(gapped.UniqueGappedKmersJob):

    _name = 'UniqueJunctionKmersJob'
    _category = 'preprocessing'
    _previous_job = ExtractJunctionKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueJunctionKmersJob(**kwargs)
        job.execute()

    def load_inputs(self):
        c = config.Configuration()
        self.gapped_kmers = {} 
        tracks = self.load_previous_job_results()
        for track in tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                gapped_kmers = json.load(json_file)
                for kmer in gapped_kmers:
                    k = canonicalize(kmer)
                    k = k[:c.hsize] + k[-c.hsize:]
                    if not k in self.gapped_kmers:
                        self.gapped_kmers[k] = gapped_kmers[kmer]
                        self.gapped_kmers[k]['tracks'] = {}
                    if not track in self.gapped_kmers[k]:
                        self.gapped_kmers[k]['tracks'][track] = 0
                        self.gapped_kmers[k]['loci'].update(gapped_kmers[kmer]['loci'])
                    self.gapped_kmers[k]['tracks'][track] += 1
        self.round_robin(tracks)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class JunctionKmersScoringJob(gapped.UniqueGappedKmersScoringJob):

    _name = 'JunctionKmersScoringJob'
    _category = 'preprocessing'
    _previous_job = UniqueJunctionKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = JunctionKmersScoringJob(**kwargs)
        job.execute()

    def transform(self, sequence, chrom):
        c = config.Configuration()
        t = time.time()
        l = len(sequence)
        for index in range(0, l - 100):
            if index % 100000 == 1:
                s = time.time()
                p = index / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
            half_mer = sequence[index: index + c.hsize]
            if not half_mer in self.half_mers:
                continue
            for k in range(32, 32 + 1):
                kmer = canonicalize(sequence[index: index + k])
                kmer = kmer[:c.hsize] + kmer[-c.hsize:]
                if kmer in self.kmers:
                    self.kmers[kmer]['count'] += 1
                    self.kmers[kmer]['loci'][chrom + '_' + str(index)] = {
                            'masks': {
                                'left': sequence[index - c.ksize: index],
                                'right': sequence[index + k: index + k + c.ksize]
                            }
                        }

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class FilterJunctionKmersJob(reduction.FilterLociIndicatorKmersJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'FilterJunctionKmersJob'
    _category = 'preprocessing'
    _previous_job = JunctionKmersScoringJob 

    @staticmethod
    def launch(**kwargs):
        job = FilterJunctionKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.tracks = self.load_previous_job_results()
        self.round_robin(self.tracks)
        print('filtering', len(self.tracks), 'tracks')

    def transform(self, track, track_name):
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                if kmers[kmer]['count'] <= 3:
                    interest_kmers = {}
                    self.kmers[kmer] = {}
                    self.kmers[kmer]['loci'] = kmers[kmer]['loci']
                    self.kmers[kmer]['total'] = 0
                    self.kmers[kmer]['count'] = 0
                    self.kmers[kmer]['doubt'] = 0
                    self.kmers[kmer]['source'] = kmers[kmer]['source']
                    self.kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                    self.kmers[kmer]['reference'] = kmers[kmer]['count']
                    self.kmers[kmer]['interest_masks'] = {}
                    for locus in self.kmers[kmer]['loci']:
                        self.kmers[kmer]['loci'][locus]['masks'] = {self.kmers[kmer]['loci'][locus]['masks']['left']: True, self.kmers[kmer]['loci'][locus]['masks']['right']: True}
                    for locus in self.kmers[kmer]['loci']:
                        if 'junction_' in locus:
                            self.kmers[kmer]['interest_masks'].update(self.kmers[kmer]['loci'][locus]['masks'])
            return None

    def is_kmer_returning(self, kmer):
        c = config.Configuration()
        for track in kmer['tracks']:
            t = c.tracks[track]
            for loci in kmer['loci']:
                if not 'junction' in loci:
                    start = int(loci.split('_')[1])
                    if abs(start - t.begin) < c.ksize:
                        return True
        return False

    def output_batch(self, batch):
        remove = {}
        for kmer in self.kmers:
            for locus in self.kmers[kmer]['loci'].keys():
                l = self.get_shared_masks(self.kmers[kmer]['interest_masks'], self.kmers[kmer]['loci'][locus]['masks'])
                if l == 0:
                    self.kmers[kmer]['loci'].pop(locus, None)
            self.kmers[kmer].pop('interest_masks', None)
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        self.kmers = {}
        for batch in self.load_output():
            for kmer in batch:
                self.kmers[kmer] = batch[kmer]
        #
        print(len(self.kmers), 'kmers')
        returning_kmers = {kmer: self.kmers[kmer] for kmer in self.kmers if self.is_kmer_returning(self.kmers[kmer])}
        with open(os.path.join(self.get_current_job_directory(), 'returning.json'), 'w') as json_file:
            json.dump(returning_kmers, json_file, indent = 4)
        self.kmers = {kmer: self.kmers[kmer] for kmer in self.kmers if not self.is_kmer_returning(self.kmers[kmer])}
        print(len(self.kmers), 'kmers after filtering returning kmers')
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        print(len(self.tracks))
        #visualizer.histogram(list(map(lambda kmer: len(self.kmers[kmer]['loci']), self.kmers)), 'number_of_loci', self.get_current_job_directory(), 'number of loci', 'number of kmers')
        #print(len(self.tracks))
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountJunctionKmersJob(map_reduce.FirstGenotypingJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'CountJunctionKmersJob'
    _category = 'preprocessing'
    _previous_job = FilterJunctionKmersJob 
    _counter_mode = 0

    @staticmethod
    def launch(**kwargs):
        job = CountJunctionKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        tracks = {}
        self.kmers = {}
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                self.kmers[kmer] = {}
                self.kmers[kmer]['loci'] = kmers[kmer]['loci']
                self.kmers[kmer]['count'] = 0 
                self.kmers[kmer]['doubt'] = 0 
                self.kmers[kmer]['total'] = 0 
                self.kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                for track in self.kmers[kmer]['tracks']:
                    tracks[track] = True
        with open(os.path.join(self.get_current_job_directory(), 'pre_inner_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        print(len(self.kmers), 'kmers')
        print(len(tracks), 'events with kmers')
        self.round_robin()

    def merge_count(self, kmer, tokens):
        count = tokens[0] 
        total = tokens[1] 
        canon = canonicalize(kmer)
        self.kmers[canon]['count'] += count / 2
        self.kmers[canon]['total'] += total / 2

    def reduce(self):
        c = config.Configuration()
        self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        for track in self.tracks:
            print('exporting track', track)
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class JunctionKmersIntegerProgrammingJob(gapped.GappedKmersIntegerProgrammingJob):

    _name = 'JunctionKmersIntegerProgrammingJob'
    _category = 'preprocessing'
    _previous_job = CountJunctionKmersJob
    _kmer_type = 'junction'

    @staticmethod
    def launch(**kwargs):
        job = JunctionKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        c = config.Configuration()
        print(cyan(track_name))
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                self.lp_kmers[kmer] = {
                    'loci': kmers[kmer]['loci'],
                    'type': self._kmer_type,
                    'count': kmers[kmer]['count'],
                    'doubt': kmers[kmer]['doubt'],
                    'total': kmers[kmer]['total'],
                    'tracks': kmers[kmer]['tracks'],
                    'weight': 1.0,
                    'coverage': c.coverage,
                    'reference': len(kmers[kmer]['loci']),
                }
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({ kmer: self.lp_kmers[kmer] for kmer in kmers }, json_file, indent = 4, sort_keys = True)
        return path

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

