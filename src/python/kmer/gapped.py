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
    counter,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
)

from kmer.sv import StructuralVariation, Inversion, Deletion
from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

import numpy
import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job for finding StructuralVariation breakpoints
# Algorithm: starts with a set of structural variation events and their approximate breakpoints and tries to refine them
# considers a radius of [-50, 50] around each end of the breakpoint, and for each pair of endpoints within that radius considers
# the area as the structural variations and applies it to the reference genome to generate a set of kmers. Discards those endpoints
# whose kmers do not all appear in the base genome the event was detected in. 
# Output: Reports the remaining boundary candidates with their list of associated kmers and the count of those kmers in the
# base genome.
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractGappedKmersJob(map_reduce.Job):

    _name = 'ExtractGappedKmersJob'
    _category = 'programming'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        extract_whole_genome()
        self.tracks = self.load_tracks()
        self.round_robin(self.tracks, lambda track: re.sub(r'\s+', '_', str(track).strip()).strip(), lambda track: track.end - track.begin > 1000000)

    # These kmers ARE NOT CANONICAL
    def transform(self, track, track_name):
        print(cyan(track_name))
        c = config.Configuration()
        gapped_kmers = track.extract_boundary_gapped_kmers()
        path = os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump(gapped_kmers, json_file, sort_keys = True, indent = 4)
        return 'gapped_kmers_' + track_name + '.json'

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueGappedKmersJob(map_reduce.Job):

    _name = 'UniqueGappedKmersJob'
    _category = 'programming'
    _previous_job = ExtractGappedKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        #self.gapped_counts_provider = counttable.JellyfishCountsProvider(self.get_gapped_reference_counts_provider())
        self.gapped_kmers = {'inner': {}, 'outer': {}}
        tracks = self.load_previous_job_results()
        for track in tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                gapped_kmers = json.load(json_file)
                for side in self.gapped_kmers:
                    for kmer in gapped_kmers[side]:
                        k = canonicalize(kmer)
                        k = k[:15] + k[-15:]
                        if not k in self.gapped_kmers[side]:
                            self.gapped_kmers[side][k] = {'tracks': {}, 'indicators': gapped_kmers[side][kmer]['indicators']}
                        if track in self.gapped_kmers[side][k]:
                            self.gapped_kmers[side][k]['tracks'][track] += 1
                        else:
                            self.gapped_kmers[side][k]['tracks'][track] = 1
        self.round_robin(tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            _gapped_kmers = json.load(json_file)
        gapped_kmers = {'inner': {}, 'outer': {}}
        for side in ['inner', 'outer']:
            for kmer in _gapped_kmers[side]:
                if 'N' in kmer:
                    continue
                k = canonicalize(kmer)
                k = k[:15] + k[-15:]
                unique = 0
                unique += len(self.gapped_kmers['inner'][k]['tracks']) if k in self.gapped_kmers['inner'] else 0
                unique += len(self.gapped_kmers['outer'][k]['tracks']) if k in self.gapped_kmers['outer'] else 0
                # The kmer is not shared with a kmer of opposite side of another event
                if unique == 1:
                    gapped_kmers[side][k] = self.gapped_kmers[side][k]
        path = os.path.join(self.get_current_job_directory(), 'unique_gapped_kmers_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump(gapped_kmers, json_file, sort_keys = True, indent = 4)
        return 'unique_gapped_kmers_' + track_name  + '.json'

    def plot(self, tracks):
        x = {'inner': [], 'outer': []}
        y = {'inner': [], 'outer': []}
        print(len(tracks))
        for track in tracks:
            print(track)
            with open(os.path.join(self.get_current_job_directory(), tracks[track]), 'r') as json_file:
                gapped_kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    print(side)
                    x[side].append(len(gapped_kmers[side]))
        visualizer.histogram(x = x['inner'], name = 'number of unique inner gapped kmers', x_label = 'number of unique gapped kmers', y_label = 'number of events', path = self.get_current_job_directory())
        visualizer.histogram(x = x['outer'], name = 'number of unique outer gapped kmers', x_label = 'number of unique gapped kmers', y_label = 'number of events', path = self.get_current_job_directory())

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueGappedKmersScoringJob(map_reduce.Job):

    _name = 'UniqueGappedKmersScoringJob'
    _category = 'programming'
    _previous_job = UniqueGappedKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueGappedKmersScoringJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.chroms = extract_whole_genome()
        self.tracks = self.load_previous_job_results()
        for track in self.tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), self.tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    for kmer in kmers[side]:
                        if len(kmer) != 30:
                            print(red(kmer))
                        self.kmers[kmer] = kmers[side][kmer]
                        self.kmers[kmer]['count'] = 0
                        self.kmers[kmer]['side'] = side
        for chrom in self.chroms:
            print(sys.getsizeof(self.chroms[chrom]))
        self.round_robin(self.chroms)

    def transform(self, sequence, chrom):
        c = config.Configuration()
        t = time.time()
        for k in range(30, 41):
            index = 0
            for kmer in stream_canonical_kmers(k, sequence):
                index += 1
                kmer = kmer[:15] + kmer[-15:]
                if kmer in self.kmers:
                    if self.kmers[kmer]['side'] == 'inner' and k != 35:
                        continue
                    self.kmers[kmer]['count'] += 1
                if index % 10000 == 0:
                    s = time.time()
                    p = (len(sequence) - index) / float(len(sequence))
                    e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                    print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))

    def output_batch(self, batch):
        n = 0
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        self.kmers = self.merge_counts()
        self.tracks = self.load_previous_job_results()
        for track in self.tracks.keys():
            with open(os.path.join(self.get_previous_job_directory(), self.tracks[track]), 'r') as json_file:
                self.tracks[track] = json.load(json_file)
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                self.tracks[track][self.kmers[kmer]['side']][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'gapped_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class SelectUniqueGappedKmersJob(counter.BaseExactCountingJob):

    _name = 'SelectUniqueGappedKmersJob'
    _category = 'programming'
    _previous_job = UniqueGappedKmersScoringJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = SelectUniqueGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        tracks = self.load_previous_job_results()
        x = []
        y = []
        self.half_mers = {}
        for track in tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                n = 0
                m = 0
                for kmer in kmers['inner']:
                    if kmers['inner'][kmer]['count'] == sum(map(lambda track: kmers['inner'][kmer]['tracks'][track], kmers['inner'][kmer]['tracks'])):
                        left = kmer[:15]
                        right = kmer[-15:]
                        self.kmers[kmer] = {'tracks': kmers['inner'][kmer]['tracks'], 'side': 'inner', 'count': {}, 'masks': kmers['inner'][kmer]['indicators']}
                        if not left in self.half_mers:
                            self.half_mers[left] = {}
                        self.half_mers[left][right] = kmer 
                        left = reverse_complement(left)
                        right = reverse_complement(right)
                        if not right in self.half_mers:
                            self.half_mers[right] = {}
                        self.half_mers[right][left] = kmer
                        n += 1
                for kmer in kmers['outer']:
                    if kmers['outer'][kmer]['count'] == 0:
                        left = kmer[:15]
                        right = kmer[-15:]
                        self.kmers[kmer] = {'tracks': kmers['outer'][kmer]['tracks'], 'side': 'outer', 'count': {}, 'masks': kmers['outer'][kmer]['indicators']}
                        if not left in self.half_mers:
                            self.half_mers[left] = {}
                        self.half_mers[left][right] = kmer 
                        left = reverse_complement(left)
                        right = reverse_complement(right)
                        if not right in self.half_mers:
                            self.half_mers[right] = {}
                        self.half_mers[right][left] = kmer 
                        m += 1
                if n + m == 0:
                    print(red('no gapped kmers found for', track))
                x.append(n)
                y.append(m)
        visualizer.histogram(x = x, name = 'number of unique inner isomers', x_label = 'percentage of novel or unique isomers', y_label = 'number of kmers', path = self.get_current_job_directory())
        visualizer.histogram(x = y, name = 'number of novel outer isomers', x_label = 'percentage of novel or unique isomers', y_label = 'number of kmers', path = self.get_current_job_directory())
        self.round_robin()

    def process_read(self, read, name, index):
        c = config.Configuration()
        kmers = index_kmers(15, read)
        for kmer in kmers:
            if kmer in self.half_mers:
                for other in self.half_mers[kmer]:
                    if not other in kmers:
                        continue
                    full = self.half_mers[kmer][other]
                    for a in kmers[kmer]:
                        for b in kmers[other]:
                            d = b - (a + 15)
                            if d >= 0:
                                if d >= 0: #list(filter(lambda mask: is_canonical_subsequence(mask, read[:a][-25:]) or is_canonical_subsequence(mask, read[b + 15:][:25]), self.kmers[full]['masks'])):
                                    if self.kmers[full]['side'] == 'outer':
                                        if d >= 0 and d <= 10:
                                            if not d in self.kmers[full]['count']:
                                                self.kmers[full]['count'][d] = 0
                                            self.kmers[full]['count'][d] += 1
                                    else:
                                        if d == 5:
                                            if not d in self.kmers[full]['count']:
                                                self.kmers[full]['count'][d] = 0
                                            self.kmers[full]['count'][d] += 1

    def reduce(self):
        c = config.Configuration()
        self.kmers = {}
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    if kmer in self.kmers:
                        for count in batch[kmer]['count']:
                            if not count in self.kmers[kmer]['count']:
                                self.kmers[kmer]['count'][count] = 0
                            self.kmers[kmer]['count'][count] += batch[kmer]['count'][count]
                    else:
                        self.kmers[kmer] = batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        # output kmers per track
        self.tracks = {}
        print(len(self.kmers))
        keys = list(self.kmers.keys())
        for kmer in keys:
            side = self.kmers[kmer]['side']
            if not self.kmers[kmer]['count']:
                self.kmers[kmer]['gap'] = -1
                self.kmers[kmer]['count'] = 0
                self.kmers[kmer]['total'] = 0
                self.kmers[kmer]['actual_gap'] = -1
            else:
                gap = max(self.kmers[kmer]['count'].items(), key = operator.itemgetter(1))[0]
                total = sum(map(lambda x: self.kmers[kmer]['count'][x], self.kmers[kmer]['count']))
                self.kmers[kmer]['gap'] = gap
                self.kmers[kmer]['count'] = self.kmers[kmer]['count'][gap]
                self.kmers[kmer]['total'] = total
                self.kmers[kmer]['actual_gap'] = -1
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner': {}, 'outer': {}}
                self.tracks[track][side][kmer] = self.kmers[kmer]
        #####################################################################################
        #if c.simulation:
        #    self.homozygous = sorted(bed.load_tracks_from_file(os.path.join(self.get_simulation_directory(), 'homozygous.bed'), ['left_drift', 'right_drift']), key = lambda track: track.begin)
        #    self.heterozygous = sorted(bed.load_tracks_from_file(os.path.join(self.get_simulation_directory(), 'heterozygous.bed'), ['left_drift', 'right_drift']), key = lambda track: track.begin)
        #    n = 0
        #    x = []
        #    wrong = 0
        #    for t in self.heterozygous:
        #        track = str(t)
        #        if not track in self.tracks:
        #            continue
        #        print(track)
        #        for kmer in self.tracks[track]['outer']:
        #            n += 1
        #            gap = str(2 + t.left_drift + 3 - t.right_drift)
        #            self.tracks[track]['outer'][kmer]['actual_gap'] = gap
        #            if self.tracks[track]['outer'][kmer]['gap'] != gap:
        #                wrong += 1
        #    for t in self.homozygous:
        #        track = str(t)
        #        if not track in self.tracks:
        #            continue
        #        print(track)
        #        for kmer in self.tracks[track]['outer']:
        #            n += 1
        #            gap = str(2 + t.left_drift + 3 - t.right_drift)
        #            self.tracks[track]['outer'][kmer]['actual_gap'] = gap
        #            if self.tracks[track]['outer'][kmer]['gap'] != gap:
        #                wrong += 1
        #    print(n, wrong)
        #####################################################################################
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'gapped_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountUniqueGappedKmersJob(map_reduce.FirstGenotypingJob, SelectUniqueGappedKmersJob):

    _name = 'CountUniqueGappedKmersJob'
    _category = 'programming'
    _previous_job = SelectUniqueGappedKmersJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountUniqueGappedKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        tracks = self.load_previous_job_results()
        self.half_mers = {}
        bed_file_name = c.bed_file.split('/')[-1]
        n = 0
        for track in tracks:
            n += 1
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    for kmer in kmers[side]:
                        # Filter bad outer kmers
                        if side == 'outer':
                            if int(kmers[side][kmer]['count']) < 10:
                                continue
                        if kmers[side][kmer]['gap'] != -1:
                            left = kmer[:15]
                            right = kmer[-15:]
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
        with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
            json.dump(self.half_mers, json_file, indent = 4)
        print(len(self.kmers), 'kmers')
        with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        self.round_robin()

    def transform(self):
        c = config.Configuration()
        if c.accelerate:
            output = subprocess.call(os.path.join(os.getcwd(), "kmer", "c_counter.out") + " " + str(self.index) + " " + self.get_current_job_directory() +  " " + c.fastq_file + " " + str(self.num_threads) + " " + "1", shell = True)
            exit()
        else:
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
                self.kmers = json.load(json_file)
            for read, name, index in self.parse_fastq():
                self.process_read(read, name, index)

    def process_read(self, read, name, index):
        c = config.Configuration()
        kmers = index_kmers(15, read)
        for kmer in kmers:
            if kmer in self.half_mers:
                for other in self.half_mers[kmer]:
                    if not other in kmers:
                        continue
                    full = self.half_mers[kmer][other]
                    for a in kmers[kmer]:
                        for b in kmers[other]:
                            d = b - (a + 15)
                            if str(d) == self.kmers[full]['gap']:
                                self.kmers[full]['count'] += 1
                                prefix = read[:a][-25:]
                                suffix = read[b + 15:][:25]
                                debug_print(kmer)
                                debug_print(other)
                                debug_print(full)
                                #debug_print(prefix)
                                #debug_print(suffix)
                                debug_print(index)
                                debug_breakpoint()
                                #if list(filter(lambda mask: is_canonical_subsequence(mask, read[:a][-25:]) or is_canonical_subsequence(mask, read[b + 15:][:25]), self.kmers[full]['masks'])):
                                #    self.kmers[full]['doubt'] += 1
                                #else:
                                #    #print(self.kmers[full]['masks'])
                                #    #print([reverse_complement(x) for x in self.kmers[full]['masks']])
                                #    #print(green(prefix), kmer, other, green(suffix))
                                #    #print(read)
                                #    debug_breakpoint()

    def reduce(self):
        c = config.Configuration()
        self.kmers = {}
        if c.accelerate:
            for i in range(0, self.num_threads):
                print('adding batch', i)
                path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
                with open (path, 'r') as json_file:
                    line = json_file.readline()
                    while line:
                        kmer = line[:line.find(':')]
                        count = int(line[line.find(':') + 1:])
                        canon = canonicalize(kmer)
                        self.kmers[canon]['count'] += count / 2
                        line = json_file.readline()
            print('exporting C kmers...')
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        else:
            for i in range(0, self.num_threads):
                print('adding batch', i)
                path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
                with open (path, 'r') as json_file:
                    batch = json.load(json_file)
                    for kmer in batch:
                        if kmer in self.kmers:
                            self.kmers[kmer]['count'] += batch[kmer]['count']
                            self.kmers[kmer]['doubt'] += batch[kmer]['doubt']
                        else:
                            self.kmers[kmer] = batch[kmer]
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        # output kmers per track
        self.tracks = {}
        for kmer in self.kmers:
            side = self.kmers[kmer]['side']
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner': {}, 'outer': {}}
                self.tracks[track][side][kmer] = self.kmers[kmer]
        #####################################################################################
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'gapped_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class GappedKmersIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'GappedKmersIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = CountUniqueGappedKmersJob
    _kmer_type = 'gapped'

    @staticmethod
    def launch(**kwargs):
        job = GappedKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            n = len(kmers['inner']) + len(kmers['outer'])
            for side in ['outer']:
                for kmer in kmers[side]:
                    if not kmer in self.lp_kmers:
                        self.lp_kmers[kmer] = {
                            'gap': kmers[side][kmer]['gap'],
                            'side': side,
                            'type': self._kmer_type,
                            'count': kmers[side][kmer]['count'],
                            'tracks': {},
                            'reference': kmers[side][kmer]['tracks'][track_name],
                            'actual_gap': kmers[side][kmer]['actual_gap'] if 'actual_gap' in kmers[side][kmer] else -1,
                        }
                    self.lp_kmers[kmer]['tracks'][track_name] = 1#kmers[side][kmer]['track']
        path = os.path.join(self.get_current_job_directory(), 'gapped_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner': {kmer: self.lp_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.lp_kmers, kmers['inner']))},
                    'outer': {kmer: self.lp_kmers[kmer] for kmer in list(filter(lambda kmer: kmer in self.lp_kmers, kmers['outer']))},
                }, json_file, indent = 4, sort_keys = True)
        return 'gapped_kmers_' + track_name + '.json'

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            kmer['residue'] = 0
            kmer['coverage'] = 42

    def generate_linear_program(self):
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
            lb = [(kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in self.lp_kmers],
            #obj = [1.0] * len(self.lp_kmers),
        )
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [1.0] * len(self.lp_kmers),
        )
        # constraints
        n = 0
        start = time.time()
        for index, kmer in enumerate(self.lp_kmers):
            if kmer['side'] == 'outer':
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
                    rhs = [kmer['count'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['G']
                )
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [-kmer['count'] + sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['G']
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
                    rhs = [kmer['count']],
                    senses = ['E']
                )
                #ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                #ind.append(len(self.tracks) + index)
                #val = list(map(lambda track: -kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                #val.append(1.0)
                #problem.linear_constraints.add(
                #    lin_expr = [cplex.SparsePair(
                #        ind = ind,
                #        val = val,
                #    )],
                #    rhs = [-kmer['count']],
                #    senses = ['G']
                #)
            self.add_error_absolute_value_constraints(problem, index)
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(self.lp_kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
        return problem

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
