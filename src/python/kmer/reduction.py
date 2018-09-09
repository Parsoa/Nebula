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
from kmer.chromosomes import *
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

class ExtractLociIndicatorKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'ExtractLociIndicatorKmersJob'
    _category = 'programming'
    _previous_job = programming.ExtractInnerKmersJob

    @staticmethod
    def launch(**kwargs):
        job = ExtractLociIndicatorKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.load_reference_counts_provider()
        self.chroms = extract_whole_genome()
        self.tracks = self.load_previous_job_results()
        self.inner_kmers = {}
        for track in self.tracks:
            print(track)
            with open(self.tracks[track], 'r') as json_file:
                kmers = json.load(json_file)['inner_kmers']
            for kmer in kmers:
                if not kmer in self.inner_kmers:
                    self.inner_kmers[kmer] = {
                        'reference': kmers[kmer]['reference'],
                        'tracks': {},
                        'loci': {}
                    }
                self.inner_kmers[kmer]['tracks'][track] = kmers[kmer]['track']
        print('Finding loci for', green(len(self.inner_kmers)), 'kmers')
        self.round_robin(self.chroms)

    def transform(self, sequence, chrom):
        c = config.Configuration()
        index = 0
        slack = (c.read_length - c.ksize) / 2
        print(cyan(chrom, len(sequence)))
        t = time.time()
        for kmer in stream_canonical_kmers(c.ksize, sequence):
            if kmer in self.inner_kmers:
                locus = chrom + '_' + str(index)
                self.inner_kmers[kmer]['loci'][locus] = {
                        'seq': {
                            'all': sequence[index - slack: index + c.ksize + slack],
                            'left': sequence[index - slack : index],
                            'right': sequence[index + c.ksize: index + c.ksize + slack]
                        }
                    }
                self.inner_kmers[kmer]['loci'][locus]['kmers'] = {
                    'left': extract_canonical_kmers(c.ksize, self.inner_kmers[kmer]['loci'][locus]['seq']['left']),
                    'right': extract_canonical_kmers(c.ksize, self.inner_kmers[kmer]['loci'][locus]['seq']['right'])
                }
            index += 1
            if index % 10000 == 0:
                s = time.time()
                p = (len(sequence) - index) / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
        return None

    def output_batch(self, batch):
        n = 0
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.inner_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            if not os.path.isfile(path):
                print(red('couldn\'t find batch'), i, red('results will be unreliable'))
                continue
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    self.inner_kmers[kmer]['loci'].update(batch[kmer]['loci'])
        for track in self.tracks:
            self.tracks[track] = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                self.tracks[track][kmer] = self.inner_kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'w') as track_file:
                json.dump(self.tracks[track], track_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'indicator_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

    def plot_kmer_reference_count(self):
        x = []
        for track in self.tracks:
            print(track)
            with open(os.path.join(self.get_current_job_directory(), self.tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)['inner_kmers']
                for kmer in kmers:
                    x.append(kmers[kmer]['reference'])
        visualizer.histogram(x = x, name = 'reference_count', path = self.get_current_job_directory(), x_label = 'kmer count in reference', y_label = 'number of kmers', step = 0.1)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountLociIndicatorKmersJob(counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'CountLociIndicatorKmersJob'
    _category = 'programming'
    _previous_job = ExtractLociIndicatorKmersJob

    @staticmethod
    def launch(**kwargs):
        job = CountLociIndicatorKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        tracks = self.load_previous_job_results()
        for track in tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                t = bed.track_from_name(track)
                for kmer in kmers:
                    seed = sum(list(map(lambda s: ord(s), kmer)))
                    if not kmer in self.kmers:
                        interest_kmers = {}
                        self.kmers[kmer] = {}
                        self.kmers[kmer]['total'] = 0
                        self.kmers[kmer]['count'] = 0
                        self.kmers[kmer]['doubt'] = 0
                        self.kmers[kmer]['tracks'] = {}
                        self.kmers[kmer]['reference'] = kmers[kmer]['reference']
                        self.kmers[kmer]['interest_kmers'] = {}
                        self.kmers[kmer]['loci'] = kmers[kmer]['loci']
                        self.kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                        for locus in self.kmers[kmer]['loci']:
                            self.kmers[kmer]['loci'][locus]['kmers'] = {k: True for k in self.kmers[kmer]['loci'][locus]['kmers']['left'].keys() + self.kmers[kmer]['loci'][locus]['kmers']['right'].keys()}
                            self.kmers[kmer]['loci'][locus]['mask'] = [
                                self.generate_kmer_mask(kmer, self.kmers[kmer]['loci'][locus]['seq']['left'], seed),
                                self.generate_kmer_mask(kmer, self.kmers[kmer]['loci'][locus]['seq']['right'], seed)
                            ]
                    n = 0
                    for locus in self.kmers[kmer]['loci']:
                        tokens = locus.split('_')
                        if tokens[0] == t.chrom and int(tokens[1]) >= t.begin and int(tokens[1]) < t.end:
                            self.kmers[kmer]['interest_kmers'].update(self.kmers[kmer]['loci'][locus]['kmers'])
                            n += 1
                    if n == 0:
                        print(yellow(track, locus))
        x = []
        for kmer in self.kmers:
            self.kmers[kmer]['indicator_kmers'] = {}
            self.kmers[kmer]['non_indicator_kmers'] = {}
            n = 0
            for locus in self.kmers[kmer]['loci'].keys():
                l = len(list(filter(lambda x: x in self.kmers[kmer]['interest_kmers'], self.kmers[kmer]['loci'][locus]['kmers'])))
                if l != 0:
                    for ukmer in self.kmers[kmer]['loci'][locus]['kmers']:
                        self.kmers[kmer]['indicator_kmers'][ukmer] = True
                    n += 1
                else:
                    for ukmer in self.kmers[kmer]['loci'][locus]['kmers']:
                        self.kmers[kmer]['non_indicator_kmers'][ukmer] = True
                    self.kmers[kmer]['loci'].pop(locus, None)
            x.append(n)
            self.kmers[kmer].pop('interest_kmers', None)
        visualizer.histogram(x = x, name = 'kmer_frequency_after_reduction', path = self.get_current_job_directory(), x_label = 'frequency', y_label = 'number of kmers', step = 1)
        print('Counting', green(len(self.kmers)), 'kmers')
        self.round_robin()

    def generate_kmer_mask(self, kmer, seq, seed):
        c = config.Configuration()
        random.seed(seed)
        indices = []
        while len(indices) != c.ksize - 2:
            i = random.randint(0, len(seq) - 1)
            if not i in indices:
                indices.append(i)
        indices = sorted(indices)
        indices.insert(0, '')
        #print(indices)
        mask = reduce(lambda x, y: x + seq[y], indices)
        #print(seq, mask)
        return mask

    # So either a locus has unique indicator kmers or it doesn't
    # if it doesn't then it means some or all of it's kmers are shared with other loci of the same kmer
    # A target either has unique kmers so the first loop should cover it or it doesn't so the second loop should
    # Read errors and SNPs when a target has unique kmers will also be handled by the second loop
    def process_read(self, read, name):
        c = config.Configuration()
        kmers = extract_canonical_kmers(c.ksize, read)
        found = False
        for kmer in kmers:
            if kmer in self.kmers:
                self.kmers[kmer]['total'] += kmers[kmer]
                for ukmer in kmers:
                    if ukmer in self.kmers[kmer]['indicator_kmers']:
                        self.kmers[kmer]['count'] += kmers[kmer]
                        found = True
                        break
                if found:
                    continue
                i = read.find(kmer)
                if i == -1:
                    i = read.find(reverse_complement(kmer))
                left = read[:i]
                right = read[i + c.ksize:]
                for locus in self.kmers[kmer]['loci']:
                    for mask in self.kmers[kmer]['loci'][locus]['mask']:
                        if is_subsequence(mask, left) or is_subsequence(mask, right):
                            self.kmers[kmer]['doubt'] += kmers[kmer]
                            break

    def reduce(self):
        if not self.resume_from_reduce:
            self.kmers = self.merge_counts('total', 'doubt')
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        else:
            self.kmers = json.load(open(os.path.join(self.get_current_job_directory(), 'kmers.json')))
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'indicator_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)
        ##############################################################################################
        x = []
        y = []
        for kmer in self.kmers:
            t = len(self.kmers[kmer]['loci'])
            x.append(self.kmers[kmer]['count'] / t)
            y.append(t)
        visualizer.histogram(x = x, name = 'average coverage after reduction', path = self.get_current_job_directory(), x_label = 'average coverage', y_label = 'number of kmers')
        visualizer.histogram(x = y, name = 'loci of interest', path = self.get_current_job_directory(), x_label = 'number of loci', y_label = 'number of kmers')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class LociIndicatorKmersIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'LociIndicatorKmersIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = CountLociIndicatorKmersJob
    _kmer_type = ''

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = LociIndicatorKmersIntegerProgrammingJob(kmer_type = 'non_unique_inner', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        print(green(track_name))
        t = bed.track_from_name(track_name)
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers.keys():
                k = self.counts_provider.get_kmer(kmer)
                if k and track_name in k['tracks']:
                    if kmers[kmer]['reference'] > 10:
                        kmers.pop(kmer, None)
                        continue
                    n = 0
                    for locus in k['loci']:
                        tokens = locus.split('_')
                        if tokens[0] == t.chrom and int(tokens[1]) >= t.begin and int(tokens[1]) < t.end:
                            n += 1
                    self.lp_kmers[kmer] = {}
                    self.lp_kmers[kmer]['type'] = self.kmer_type
                    self.lp_kmers[kmer]['count'] = k['count']
                    self.lp_kmers[kmer]['doubt'] = k['doubt']
                    self.lp_kmers[kmer]['total'] = k['total']
                    self.lp_kmers[kmer]['reference'] = len(kmers[kmer]['loci'])
                    self.lp_kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                else:
                    kmers.pop(kmer, None)
            if not kmers: 
                print('no inner kmers found for', red(track_name))
                return None
        path = os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in kmers}, json_file, indent = 4, sort_keys = True)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class DebugCountLociIndicatorKmersJob(counter.SimulationExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = DebugCountLociIndicatorKmersJob(job_name = 'DebugCountLociIndicatorKmersJob_', previous_job_name = 'CountLociIndicatorKmersJob_', category = 'programming', batch_file_prefix = 'indicator_kmers', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
        self.homozygous = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'homozygous.bed'))], key = lambda track: track.begin)
        self.homozygous_t = {re.sub(r'\s+', '_', str(track).strip()).strip(): True for track in self.homozygous}
        self.heterozygous = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'heterozygous.bed'))], key = lambda track: track.begin)
        self.heterozygous_t = {re.sub(r'\s+', '_', str(track).strip()).strip(): True for track in self.heterozygous} 
        self.present = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'present.bed'))], key = lambda track: track.begin)
        self.kmers = {}
        for kmer in kmers:
            if len(kmers[kmer]['tracks']) != 1:
                continue
            self.kmers[kmer] = {}
            self.kmers[kmer]['count'] = kmers[kmer]['count']
            self.kmers[kmer]['total'] = 0
            self.kmers[kmer]['doubt'] = 0
            self.kmers[kmer]['reference'] = kmers[kmer]['reference']
            self.kmers[kmer]['loci'] = {'1': {}, '2': {}}
            self.kmers[kmer]['actual_loci'] = []
            for track in kmers[kmer]['tracks']:
                t = bed.track_from_name(track)
                for locus in kmers[kmer]['loci'].keys():
                    tokens = locus.split('_')
                    if tokens[0] == t.chrom and int(tokens[1]) >= t.begin and int(tokens[1]) <= t.end - c.ksize:
                        locus_1, locus_2 = self.convert_loci(locus, track)
                        if locus_1:
                            self.kmers[kmer]['loci']['1'][locus_1] = 0
                        if locus_2:
                            self.kmers[kmer]['loci']['2'][locus_2] = 0
                        self.kmers[kmer]['actual_loci'].append(locus)
                    else:
                        self.kmers[kmer]['loci'].pop(locus, None)
                self.kmers[kmer]['track'] = track
                self.kmers[kmer]['zygosity'] = 'homozygous' if track in self.homozygous_t\
                    else 'heterozygous' if track in self.heterozygous_t\
                    else 'absent'
        self.round_robin()

    def convert_loci(self, locus, track):
        tokens = locus.split('_')
        if tokens[0] != 'chr2':
            return locus, locus
        begin = int(tokens[1])
        offset_1 = sum(list(map(lambda y: y.end - y.begin, list(filter(lambda x: x.end < begin, self.present)))))
        offset_2 = sum(list(map(lambda y: y.end - y.begin, list(filter(lambda x: x.end < begin, self.homozygous)))))
        if track in self.homozygous_t:
            return None, None
        if track in self.heterozygous_t:
            return None, tokens[0] + '_' + str(begin - offset_2)
        else:
            return tokens[0] + '_' + str(begin - offset_1), tokens[0] + '_' + str(begin - offset_2)

    def process_read(self, read, name, track_name):
        c = config.Configuration()
        name = name[1:name.find(':')]
        strand = '1' if '_1.1' in track_name else '2'
        tokens = name.split('_')
        chrom = tokens[0]
        begin = int(tokens[1])
        end = int(tokens[2])
        kmers = extract_canonical_kmers(c.ksize, read)
        for kmer in kmers:
            if kmer in self.kmers:
                found = False
                for locus in self.kmers[kmer]['loci'][strand]:
                    tokens = locus.split('_')
                    if tokens[0] == chrom:
                        if int(tokens[1]) >= begin and int(tokens[1]) < end:
                            self.kmers[kmer]['loci'][strand][locus] += kmers[kmer]
                            break
                self.kmers[kmer]['total'] += kmers[kmer]

    def reduce(self):
        self.tracks = {}
        # Begin temporary fix
        for kmer in self.kmers.keys():
            if not 'track' in self.kmers[kmer]:
                self.kmers.pop(kmer, None)
        # End temporary fix
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    # Begin temporary fix
                    if not kmer in self.kmers:
                        continue
                    # End temporary fix
                    self.kmers[kmer]['total'] += batch[kmer]['total']
                    for strand in self.kmers[kmer]['loci']:
                        for locus in self.kmers[kmer]['loci'][strand]:
                            self.kmers[kmer]['loci'][strand][locus] += batch[kmer]['loci'][strand][locus]
                    track = self.kmers[kmer]['track']
                    if not track in self.tracks:
                        self.tracks[track] = {}
                    self.tracks[track][kmer] = self.kmers[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_merge.json'), 'w') as json_file:
            json.dump({track: self.batch_file_prefix + '_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class DebugLociIndicatorKmersIntegerProgrammingJob(programming.IntegerProgrammingJob):

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = DebugLociIndicatorKmersIntegerProgrammingJob(job_name = 'DebugLociIndicatorKmersIntegerProgrammingJob_', previous_job_name = 'DebugCountLociIndicatorKmersJob_', category = 'programming', batch_file_prefix = 'indicator_kmers', previous_job_batch_file_prefix = 'indicator_kmers', kmer_type = 'non_unique_inner', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        t = bed.track_from_name(track_name)
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers.keys():
                k = self.counts_provider.get_kmer(kmer)
                if k and k['track'] == track_name:
                    self.lp_kmers[kmer] = {}
                    self.lp_kmers[kmer]['count'] = sum(map(lambda side: sum(map(lambda locus: k['loci'][side][locus], k['loci'][side])), k['loci']))
                    self.lp_kmers[kmer]['total'] = k['total']
                    self.lp_kmers[kmer]['reference'] = len(kmers[kmer]['actual_loci'])
                    self.lp_kmers[kmer]['tracks'] = {track_name: len(kmers[kmer]['actual_loci'])}
                else:
                    kmers.pop(kmer, None)
            if not kmers: 
                print('no inner kmers found for', red(track_name))
                return None
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    self.kmer_type + '_kmers': {kmer: self.lp_kmers[kmer] for kmer in kmers},
                }, json_file, indent = 4, sort_keys = True)
        return path 

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
