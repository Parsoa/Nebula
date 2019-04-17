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
    _category = 'programming'
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

    def load_inputs(self):
        c = config.Configuration()
        self.bamfile = pysam.AlignmentFile(c.bam, "rb")
        self.tracks = c.tracks
        self.extract_kmers()
        self.round_robin({})
        self.num_threads = 1

    def extract_kmers(self):
        c = config.Configuration()
        contig_index = {}
        output = {}
        for track in self.tracks:
            contig_index[self.tracks[track].contig] = self.tracks[track]
        n = 0
        print(len(contig_index), 'total tracks')
        for read in self.bamfile.fetch():
            t = contig_index[read.query_name]
            path = self.transform(t, read.query_sequence)
            if path:
                output[t.id] = path
                n += 1
            if n % 1000 == 0:
                print(n)
        with open(os.path.join(self.get_current_job_directory(), 'no_contig.bed'), 'w') as bed_file:
            for contig in contig_index:
                track = contig_index[contig]
                bed_file.write(track.bedify())

    def transform(self, track, assembly):
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
                    'loci': {
                        'junction_' + str(i) + '_' + track.id: {
                            'masks': {
                                'left': assembly[i - c.ksize: i],
                                'right': assembly[i + c.ksize: i + c.ksize + c.ksize]
                            }
                        }
                    }
                }
            junction_kmers[kmer]['count'] += 1
        junction_kmers = {kmer: junction_kmers[kmer] for kmer in junction_kmers if junction_kmers[kmer]['count'] == 1}
        path = os.path.join(self.get_current_job_directory(), track.name + '.json')
        if len(junction_kmers) > 0:
            with open(path, 'w') as json_file:
                json.dump(junction_kmers, json_file, indent = 4)
            return track.name + '.json'
        return None

    #def transform(self, track, track_name):
    #    c = config.Configuration()
    #    junction_kmers = {}
    #    t = c.tracks[track_name]
    #    #reads = self.bamfile.fetch(contig = t.contig)
    #    command = 'samtools view ' + c.bam + ' | grep "' + t.contig + '"'
    #    output = subprocess.check_output(command, shell = True)
    #    assembly = output.split()[9]
    #    if t.svtype == 'DEL':
    #        l = list(range(int(t.contig_start) - 24, int(t.contig_start) - 8))
    #    if t.svtype == 'INS':
    #        l = list(range(int(t.contig_start) - 24, int(t.contig_start) - 8)) + list(range(int(t.contig_end) - 24, int(t.contig_end) - 8))
    #    for i in l:
    #        kmer = canonicalize(assembly[i: i + c.ksize])
    #        if not kmer in junction_kmers:
    #            junction_kmers[kmer] = {
    #                'count': 0,
    #                'loci': {
    #                    'junction_' + str(i) + '_' + track_name: {
    #                        'masks': {
    #                            'left': assembly[i - c.ksize: i],
    #                            'right': assembly[i + c.ksize: i + c.ksize + c.ksize]
    #                        }
    #                    }
    #                }
    #            }
    #        junction_kmers[kmer]['count'] += 1
    #    junction_kmers = {kmer: junction_kmers[kmer] for kmer in junction_kmers if junction_kmers[kmer]['count'] == 1}
    #    path = os.path.join(self.get_current_job_directory(), track_name + '.json')
    #    if len(junction_kmers) > 0:
    #        with open(path, 'w') as json_file:
    #            json.dump(junction_kmers, json_file, indent = 4)
    #        return track_name + '.json'
    #    return None

    #def transform(self, track, track_name):
    #    c = config.Configuration()
    #    junction_kmers = {}
    #    slack = 70
    #    stride = c.ksize
    #    self.sv_type = self.get_sv_type()
    #    if not c.simulation:
    #        track = track.lift()
    #    reads = chain(self.bamfile.fetch(track.chrom, track.begin - slack, track.begin + slack), self.bamfile.fetch(track.chrom, track.end - slack, track.end + slack))
    #    n = 0
    #    strid = c.ksize
    #    for read in reads:
    #        n += 1
    #        if read.query_alignment_length != len(read.query_sequence):
    #            seq = read.query_sequence
    #            if read.reference_start >= track.begin - slack and read.reference_start <= track.end + slack:
    #                cigar = read.cigartuples
    #                clips = []
    #                offset = 0
    #                for c in cigar:
    #                    if c[0] == 4: #soft clip
    #                        clips.append((offset, offset + c[1]))
    #                    offset += c[1]
    #                index = 0
    #                for kmer in stream_kmers(stride, False, True, seq):
    #                    if self.is_clipped((index, index + stride), clips):
    #                        if not kmer in junction_kmers:
    #                            junction_kmers[kmer] = {
    #                                'count': 0,
    #                                'loci': {
    #                                    'break_point_' + track_name: {
    #                                        'masks': {
    #                                            'left': None,
    #                                            'right': None
    #                                        }
    #                                    }
    #                                }
    #                            }
    #                        if len(seq[index - stride: index]) == stride:
    #                            junction_kmers[kmer]['loci']['break_point_' + track_name]['masks']['left'] = seq[index - stride: index]
    #                        if len(seq[index + stride: index + stride + stride]) == stride:
    #                            junction_kmers[kmer]['loci']['break_point_' + track_name]['masks']['right'] = seq[index + stride: index + stride + stride]
    #                        junction_kmers[kmer]['count'] += 1
    #                    index += 1
    #    _junction_kmers = {} 
    #    for kmer in junction_kmers:
    #        if junction_kmers[kmer]['count'] >= 3 and junction_kmers[kmer]['loci']['break_point_' + track_name]['masks']['right'] and junction_kmers[kmer]['loci']['break_point_' + track_name]['masks']['left']:
    #            _junction_kmers[kmer] = junction_kmers[kmer]
    #    path = os.path.join(self.get_current_job_directory(), track_name + '.json')
    #    if len(_junction_kmers) > 0:
    #        with open(path, 'w') as json_file:
    #            json.dump(_junction_kmers, json_file, indent = 4)
    #        return track_name + '.json'
    #    return None

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
    _category = 'programming'
    _previous_job = ExtractJunctionKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueJunctionKmersJob(**kwargs)
        job.execute()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class JunctionKmersScoringJob(gapped.UniqueGappedKmersScoringJob):

    _name = 'JunctionKmersScoringJob'
    _category = 'programming'
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
            for k in range(32, 33):
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
    _category = 'programming'
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

    def transform(self, track, track_name):
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            t = bed.track_from_name(track_name)
            for kmer in kmers:
                if kmers[kmer]['count'] <= 3:
                    interest_kmers = {}
                    self.kmers[kmer] = {}
                    self.kmers[kmer]['loci'] = kmers[kmer]['loci']
                    self.kmers[kmer]['total'] = 0
                    self.kmers[kmer]['count'] = 0
                    self.kmers[kmer]['doubt'] = 0
                    self.kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                    self.kmers[kmer]['reference'] = kmers[kmer]['count']
                    self.kmers[kmer]['interest_masks'] = {}
                    for locus in self.kmers[kmer]['loci']:
                        self.kmers[kmer]['loci'][locus]['masks'] = {self.kmers[kmer]['loci'][locus]['masks']['left']: True, self.kmers[kmer]['loci'][locus]['masks']['right']: True}
                    for locus in self.kmers[kmer]['loci']:
                        if 'break_point_' in locus:
                            self.kmers[kmer]['interest_masks'].update(self.kmers[kmer]['loci'][locus]['masks'])
            return None

    def output_batch(self, batch):
        for kmer in self.kmers:
            for locus in self.kmers[kmer]['loci'].keys():
                l = self.get_shared_masks(self.kmers[kmer]['interest_masks'], self.kmers[kmer]['loci'][locus]['masks'])
                if l == 0:
                    self.kmers[kmer]['loci'].pop(locus, None)
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

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
    _category = 'programming'
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
    _category = 'programming'
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

