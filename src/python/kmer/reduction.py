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

    @staticmethod
    def launch(**kwargs):
        job = ExtractLociIndicatorKmersJob(job_name = 'ExtractLociIndicatorKmersJob_', previous_job_name = '', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.bedtools = {str(track): track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'all.bed'))}
        self.round_robin(self.bedtools, lambda track: re.sub(r'\s+', '_', str(track).strip()).strip(), lambda track: track.end - track.start > 1000000)
        self.chroms = extract_whole_genome()

    def transform(self, track, track_name):
        c = config.Configuration()
        sv = self.get_sv_type()(bed.track_from_name(track_name))
        _inner_kmers = sv.get_inner_kmers(counter = self.reference_counts_provider.get_kmer_count, count = 10, n = 1000)
        novel_kmers = sv.get_boundary_kmers(begin = 0, end = 0, counter = self.reference_counts_provider.get_kmer_count, count = 0)
        l = len(_inner_kmers)
        _inner_kmers = {kmer: _inner_kmers[kmer] for kmer in filter(lambda k: k not in novel_kmers and reverse_complement(k) not in novel_kmers, _inner_kmers)}
        if l != len(_inner_kmers):
            print(yellow(track_name))
        inner_kmers = {}
        for kmer in _inner_kmers:
            inner_kmers[kmer] = {}
            inner_kmers[kmer]['loci'] = self.find_kmer_loci(kmer, track_name)
            inner_kmers[kmer]['reference'] = self.reference_counts_provider.get_kmer_count(kmer)
            inner_kmers[kmer]['track'] = _inner_kmers[kmer]
        path = os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track_name  + '.json') 
        with open(path, 'w') as json_file:
            json.dump({'inner_kmers': inner_kmers}, json_file, sort_keys = True, indent = 4)
        return 'inner_kmers_' + track_name  + '.json'

    def find_kmer_loci(self, kmer, track_name):
        loci = {}
        track = bed.track_from_name(track_name)
        for chrom in self.chroms:
            t = self.find_all(self.chroms[chrom], kmer)
            t += self.find_all(self.chroms[chrom], reverse_complement(kmer))
            for position in t:
                name = chrom + '_' + str(position)
                slack = (c.read_length - c.ksize) / 2
                loci[name] = {}
                loci[name]['seq'] = {}
                loci[name]['seq']['all'] = self.chroms[chrom][position - slack : position + c.ksize + slack]
                loci[name]['seq']['left'] = self.chroms[chrom][position - slack : position]
                loci[name]['seq']['right'] = self.chroms[chrom][position + c.ksize : position + c.ksize + slack]
                loci[name]['kmers'] = {}
                loci[name]['kmers']['left'] = extract_canonical_kmers(loci[name]['seq']['left'])
                loci[name]['kmers']['right'] = extract_canonical_kmers(loci[name]['seq']['right'])
                #o[name]['unique_kmers'] = { 'right': {}, 'left': {} }
        #kmers = {}
        #for i in o:
        #    for side in o[i]['kmers']:
        #        for kmer in o[i]['kmers'][side]:
        #            if not kmer in kmers:
        #                kmers[kmer] = {}
        #            kmers[kmer][i] = True
        #for i in o:
        #    for side in o[i]['kmers']:
        #        for kmer in o[i]['kmers'][side]:
        #            if kmer in kmers and len(kmers[kmer]) == 1:
        #                o[i]['unique_kmers'][side][kmer] = True
        return loci

    def find_all(self, string, substring):
        l = []
        index = -1
        while True:
            index = string.find(substring, index + 1)
            if index == -1:  
                break
            l.append(index)
        return l

    #def reduce(self):
    #    tracks = map_reduce.Job.reduce(self)
    #    x = []
    #    y = []
    #    z = []
    #    bad = {}
    #    with open(os.path.join('/share/hormozdiarilab/Codes/NebulousSerendipity/simulation/HG00513.hg38.chr2.DEL.bed/31/MixGappedInnerKmersIntegerProgrammingJob', 'no_signal.bed')) as bed_file:
    #        lines = bed_file.readlines()
    #        for line in lines:
    #            tokens = line.split()
    #            name = tokens[0] + '_' + tokens[1] + '_' + tokens[2]
    #            bad[name] = True
    #    with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
    #        for track in tracks:
    #            print(track)
    #            if not track in bad:
    #                continue
    #            with open(os.path.join(self.get_current_job_directory(), tracks[track]), 'r') as json_file:
    #                kmers = json.load(json_file)
    #                for kmer in kmers['inner_kmers']:
    #                    c = kmers['inner_kmers'][kmer]['count']
    #                    if c > 1:
    #                        n = len(list(filter(lambda o: len(kmers['inner_kmers'][kmer]['ocurrences'][o]['unique_kmers']['left']) > 0 or len(kmers['inner_kmers'][kmer]['ocurrences'][o]['unique_kmers']['right']) > 0, kmers['inner_kmers'][kmer]['ocurrences'])))
    #                        x.append(c - n)
    #                l = len(kmers['inner_kmers'])
    #                m = len(list(filter(lambda kmer: len(list(filter(lambda o: len(kmers['inner_kmers'][kmer]['ocurrences'][o]['unique_kmers']['left']) > 0 or len(kmers['inner_kmers'][kmer]['ocurrences'][o]['unique_kmers']['right']) > 0, kmers['inner_kmers'][kmer]['ocurrences']))) == len(kmers['inner_kmers'][kmer]['ocurrences']), kmers['inner_kmers'])))
    #                t = bed.track_from_name(track)
    #                bed_file.write(t.chrom + '\t' + str(t.start) + '\t' + str(t.end) + '\t' + str(len(kmers['unique_inner_kmers'])) + '\t' +
    #                    str(l) + '\t' +
    #                    str(m) + '\n')
    #                if l != 0:
    #                    z.append(float(m) / float(l))
    #                y.append(l)
    #    visualizer.histogram(x = x, name = 'kmer_dimension_reduction', path = self.get_current_job_directory(), x_label = 'kmer frequency after reduction', y_label = 'number of kmers', step = 1)
    #    visualizer.histogram(x = z, name = 'percentage_local_unique_kmers', path = self.get_current_job_directory(), x_label = 'percentage of inner kmers with unique markers', y_label = 'number of events', step = 0.1)
    #    visualizer.histogram(x = y, name = 'num_inner_kmers', path = self.get_current_job_directory(), x_label = 'number of inner kmers', y_label = 'number of events', step = 0.1)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountLociIndicatorKmersJob(counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountLociIndicatorKmersJob(job_name = 'CountLociIndicatorKmersJob_', previous_job_name = 'ExtractLociIndicatorKmersJob_', category = 'programming', **kwargs)
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
        a = 0
        b = 0
        for track in tracks:
            print(cyan(track))
            with open(os.path.join(self.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                t = bed.track_from_name(track)
                for kmer in kmers['inner_kmers']:
                    a += 1
                    seed = sum(list(map(lambda s: ord(s), kmer)))
                    if not kmer in self.kmers:
                        interest_kmers = {}
                        self.kmers[kmer] = {}
                        self.kmers[kmer]['total'] = 0
                        self.kmers[kmer]['count'] = 0
                        self.kmers[kmer]['doubt'] = 0
                        self.kmers[kmer]['track'] = track
                        self.kmers[kmer]['reference'] = kmers['inner_kmers'][kmer]['count']
                        self.kmers[kmer]['indicator_kmers'] = {}
                        y.append(self.kmers[kmer]['reference'])
                        self.kmers[kmer]['loci'] = kmers['inner_kmers'][kmer]['ocurrences']
                        for locus in self.kmers[kmer]['loci']:
                            self.kmers[kmer]['loci'][locus]['kmers'] = {k: True for k in self.kmers[kmer]['loci'][locus]['kmers']['left'].keys() + self.kmers[kmer]['loci'][locus]['kmers']['right'].keys()}
                            #self.kmers[kmer]['loci'][locus]['unique_kmers'] = {k: True for k in self.kmers[kmer]['loci'][locus]['unique_kmers']['left'].keys() + self.kmers[kmer]['loci'][locus]['unique_kmers']['right'].keys()}
                            self.kmers[kmer]['loci'][locus]['mask'] = [
                                self.generate_kmer_mask(kmer, self.kmers[kmer]['loci'][locus]['seq']['left'], seed),
                                self.generate_kmer_mask(kmer, self.kmers[kmer]['loci'][locus]['seq']['right'], seed)
                            ]
                            tokens = locus.split('_')
                            if tokens[0] == t.chrom and int(tokens[1]) >= t.start and int(tokens[1]) < t.end:
                                interest_kmers.update(self.kmers[kmer]['loci'][locus]['kmers'])
                        n = 0
                        for locus in self.kmers[kmer]['loci'].keys():
                            l = len(list(filter(lambda x: x in interest_kmers, self.kmers[kmer]['loci'][locus]['kmers'])))
                            if l != 0:
                                for ukmer in self.kmers[kmer]['loci'][locus]['kmers']:
                                    self.kmers[kmer]['indicator_kmers'][ukmer] = True
                                n += 1
                            else:
                                self.kmers[kmer]['loci'].pop(locus, None)
                        x.append(n)
                    else:
                        #print(kmer, self.kmers[kmer]['reference'])
                        b += 1
                        self.kmers.pop(kmer, None)
        print(green(a), 'all', green(len(self.kmers)), 'remain,', green(b), 'duplicates')
        exit()
        #visualizer.histogram(x = x, name = 'number of counted loci', path = self.get_current_job_directory(), x_label = 'number of loci', y_label = 'number of kmers', step = 1)
        #visualizer.histogram(x = y, name = 'number of reference loci', path = self.get_current_job_directory(), x_label = 'number of loci', y_label = 'number of kmers', step = 1)

    def generate_kmer_mask(self, kmer, seq, seed):
        random.seed(seed)
        indices = []
        while len(indices) != len(seq) / 2:
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
                #if found:
                #    continue
                #i = read.find(kmer)
                #if i == -1:
                #    i = read.find(reverse_complement(kmer))
                #left = read[:i]
                #right = read[i + c.ksize:]
                #for locus in self.kmers[kmer]['loci']:
                #    if found:
                #        break
                #    for mask in self.kmers[kmer]['loci'][locus]['mask']:
                #        if is_subsequence(mask, left) or is_subsequence(mask, right):
                #            self.kmers[kmer]['doubt'] += kmers[kmer]
                #            found = True
                #            break

    def reduce(self):
        #if not self.resume_from_reduce:
        self.kmers = self.merge_counts('total')
        #else:
        #    self.kmers = json.load(open(os.path.join(self.get_current_job_directory(), 'kmers.json')))
        self.tracks = {}
        print(n, '/', len(self.kmers))
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
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

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = LociIndicatorKmersIntegerProgrammingJob(job_name = 'LociIndicatorKmersIntegerProgrammingJob_', previous_job_name = 'CountLociIndicatorKmersJob_', category = 'programming', batch_file_prefix = 'indicator_kmers', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        print(green(track_name))
        t = bed.track_from_name(track_name)
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            inner_kmers = json.load(json_file)
            for kmer in inner_kmers.keys():
                k = self.counts_provider.get_kmer(kmer)
                if k and k['track'] == track_name:
                    if inner_kmers[kmer]['reference'] == 1:
                        print(kmer)
                        inner_kmers.pop(kmer, None)
                        continue
                    n = 0
                    for locus in k['loci']:
                        tokens = locus.split('_')
                        if tokens[0] == t.chrom and int(tokens[1]) >= t.start and int(tokens[1]) < t.end:
                           n += 1
                    self.inner_kmers[kmer] = {}
                    self.inner_kmers[kmer]['count'] = k['count']
                    self.inner_kmers[kmer]['doubt'] = k['doubt']
                    self.inner_kmers[kmer]['total'] = k['total']
                    self.inner_kmers[kmer]['reference'] = len(inner_kmers[kmer]['loci'])
                    self.inner_kmers[kmer]['tracks'] = {track_name: n}
                else:
                    inner_kmers.pop(kmer, None)
            if not inner_kmers: 
                print('no inner kmers found for', red(track_name))
                return None
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner_kmers': {kmer: self.inner_kmers[kmer] for kmer in inner_kmers},
                }, json_file, indent = 4, sort_keys = True)
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
        job = DebugCountLociIndicatorKmersJob(job_name = 'DebugCountLociIndicatorKmersJob_', previous_job_name = 'CountLociIndicatorKmersJob_', category = 'reduction', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
        self.homozygous = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'homozygous.bed'))], key = lambda track: track.start)
        self.homozygous_t = {re.sub(r'\s+', '_', str(track).strip()).strip(): True for track in self.homozygous}
        self.heterozygous = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'heterozygous.bed'))], key = lambda track: track.start)
        self.heterozygous_t = {re.sub(r'\s+', '_', str(track).strip()).strip(): True for track in self.heterozygous} 
        self.present = sorted([track for track in pybedtools.BedTool(os.path.join(self.get_simulation_directory(), 'present.bed'))], key = lambda track: track.start)
        self.kmers = {}
        for kmer in kmers:
            track = kmers[kmer]['track']
            t = bed.track_from_name(track)
            self.kmers[kmer] = {}
            self.kmers[kmer]['count'] = kmers[kmer]['count']
            self.kmers[kmer]['total'] = 0
            self.kmers[kmer]['doubt'] = 0
            n = 0
            for locus in kmers[kmer]['loci']:
                tokens = locus.split('_')
                if tokens[0] == t.chrom and int(tokens[1]) >= t.start and int(tokens[1]) < t.end:
                    n += 1
            self.kmers[kmer]['track'] = track
            self.kmers[kmer]['reference'] = kmers[kmer]['reference']
            self.kmers[kmer]['loci'] = {'1': {}, '2': {}}
            self.kmers[kmer]['actual_loci'] = []
            for locus in kmers[kmer]['loci']:
                locus_1, locus_2 = self.convert_loci(locus, self.kmers[kmer]['track'])
                if locus_1:
                    self.kmers[kmer]['loci']['1'][locus_1] = 0
                if locus_2:
                    self.kmers[kmer]['loci']['2'][locus_2] = 0
                self.kmers[kmer]['actual_loci'].append(locus)
            self.kmers[kmer]['zygosity'] = 'homozygous' if track in self.homozygous_t else 'heterozygous' if track in self.heterozygous_t else 'absent'
        self.round_robin()

    def convert_loci(self, locus, track):
        tokens = locus.split('_')
        if tokens[0] != 'chr2':
            return locus, locus
        start = int(tokens[1])
        offset_1 = sum(list(map(lambda y: y.end - y.start, list(filter(lambda x: x.end < start, self.present)))))
        offset_2 = sum(list(map(lambda y: y.end - y.start, list(filter(lambda x: x.end < start, self.homozygous)))))
        if track in self.homozygous_t:
            return None, None
        if track in self.heterozygous_t:
            return None, tokens[0] + '_' + str(start - offset_2)
        else:
            return tokens[0] + '_' + str(start - offset_1), tokens[0] + '_' + str(start - offset_2)

    def process_read(self, read, name):
        c = config.Configuration()
        name = name[1:name.find(':')]
        strand = '1' if '_1.1' in track_name else '2'
        tokens = name.split('_')
        chrom = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])
        kmers = extract_canonical_kmers(c.ksize, read)
        for kmer in kmers:
            if kmer in self.kmers:
                found = False
                for locus in self.kmers[kmer]['loci'][strand]:
                    tokens = locus.split('_')
                    if tokens[0] == chrom:
                        if int(tokens[1]) >= start and int(tokens[1]) < end:
                            self.kmers[kmer]['loci'][strand][locus] += kmers[kmer]
                            break
                self.kmers[kmer]['total'] += kmers[kmer]

    def reduce(self):
        self.tracks = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
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
            with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'indicator_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

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
        job = DebugLociIndicatorKmersIntegerProgrammingJob(job_name = 'DebugLociIndicatorKmersIntegerProgrammingJob_', previous_job_name = 'DebugCountLociIndicatorKmersJob_', category = 'programming', batch_file_prefix = 'indicator_kmers', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        t = bed.track_from_name(track_name)
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            inner_kmers = json.load(json_file)
            for kmer in inner_kmers.keys():
                k = self.counts_provider.get_kmer(kmer)
                if k and k['track'] == track_name:
                    n = 0
                    for locus in k['actual_loci']:
                        tokens = locus.split('_')
                        if tokens[0] == t.chrom and int(tokens[1]) >= t.start and int(tokens[1]) < t.end:
                           n += 1
                    self.inner_kmers[kmer] = {}
                    self.inner_kmers[kmer]['count'] = sum(map(lambda side: sum(map(lambda locus: k['loci'][side][locus], k['loci'][x])), k['loci']))
                    self.inner_kmers[kmer]['total'] = k['total']
                    self.inner_kmers[kmer]['reference'] = len(inner_kmers[kmer]['actual_loci'])#len(inner_kmers[kmer]['loci'])
                    self.inner_kmers[kmer]['tracks'] = {track_name: n}
                else:
                    inner_kmers.pop(kmer, None)
            if not inner_kmers: 
                print('no inner kmers found for', red(track_name))
                return None
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(
                {
                    'inner_kmers': {kmer: self.inner_kmers[kmer] for kmer in inner_kmers},
                }, json_file, indent = 4, sort_keys = True)
        return path 


# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
