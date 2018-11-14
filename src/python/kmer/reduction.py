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
import subprocess

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

import scipy

from kmer.debug import *
from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

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
        self.chroms = extract_whole_genome()
        self.tracks = self.load_previous_job_results()
        self.inner_kmers = {}
        for track in self.tracks:
            print(track)
            with open(os.path.join(self.get_previous_job_directory(), self.tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
            for kmer_type in ['non_unique_inner_kmers', 'unique_inner_kmers']:
                for kmer in kmers[kmer_type]:
                    if not kmer in self.inner_kmers:
                        self.inner_kmers[kmer] = {
                            'reference': kmers[kmer_type][kmer]['reference'],
                            'tracks': {},
                            'loci': {}
                        }
                    self.inner_kmers[kmer]['tracks'][track] = kmers[kmer_type][kmer]['track']
        print('Finding loci for', green(len(self.inner_kmers)), 'kmers')
        self.round_robin(self.chroms)

    def transform(self, sequence, chrom):
        c = config.Configuration()
        index = 0
        slack = c.ksize#(c.read_length - c.ksize) / 2
        print(cyan(chrom, len(sequence)))
        t = time.time()
        for kmer in stream_kmers(c.ksize, True, sequence):
            if kmer in self.inner_kmers:
                locus = chrom + '_' + str(index)
                self.inner_kmers[kmer]['loci'][locus] = {
                    'seq': {
                        'all': sequence[index - slack: index + c.ksize + slack],
                        'left': sequence[index - slack: index],
                        'right': sequence[index + c.ksize: index + c.ksize + slack]
                    }
                }
                self.inner_kmers[kmer]['loci'][locus]['kmers'] = {
                    'left': extract_kmers(c.ksize, True, self.inner_kmers[kmer]['loci'][locus]['seq']['left']),
                    'right': extract_kmers(c.ksize, True, self.inner_kmers[kmer]['loci'][locus]['seq']['right'])
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

class FilterLociIndicatorKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'FilterLociIndicatorKmersJob'
    _category = 'programming'
    _previous_job = ExtractLociIndicatorKmersJob

    @staticmethod
    def launch(**kwargs):
        job = FilterLociIndicatorKmersJob(**kwargs)
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
        print(cyan(track))
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            t = bed.track_from_name(track_name)
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
                    self.kmers[kmer]['interest_masks'] = {}
                    self.kmers[kmer]['loci'] = kmers[kmer]['loci']
                    self.kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                    for locus in self.kmers[kmer]['loci']:
                        self.kmers[kmer]['loci'][locus]['kmers'] = {k: True for k in self.kmers[kmer]['loci'][locus]['kmers']['left'].keys() + self.kmers[kmer]['loci'][locus]['kmers']['right'].keys()}
                        self.kmers[kmer]['loci'][locus]['masks'] = self.generate_kmer_mask(kmer, self.kmers[kmer]['loci'][locus]['seq']['left'], self.kmers[kmer]['loci'][locus]['seq']['right'], seed)
                n = 0
                for locus in self.kmers[kmer]['loci']:
                    tokens = locus.split('_')
                    if tokens[0].lower() == t.chrom.lower() and int(tokens[1]) >= t.begin and int(tokens[1]) < t.end:
                        self.kmers[kmer]['interest_kmers'].update(self.kmers[kmer]['loci'][locus]['kmers'])
                        self.kmers[kmer]['interest_masks'].update(self.kmers[kmer]['loci'][locus]['masks'])
                        n += 1
                if n == 0:
                    print(yellow(track, locus))
        return None

    def generate_kmer_mask(self, kmer, left, right, seed):
        c = config.Configuration()
        return {left: True, right: True}
        random.seed(seed)
        masks = {}
        for seq in [left, right]:
            for j in range(0, 5):
                indices = {}
                mask = 'N' * len(seq)
                while len(indices) != c.ksize - 2:
                    i = random.randint(0, len(seq) - 1)
                    if not i in indices:
                        indices[i] = True
                        mask = mask[:i] + seq[i] + mask[i + 1:]
                masks[mask] = True
        print(masks)
        return masks

    def output_batch(self, batch):
        for kmer in self.kmers:
            self.kmers[kmer]['indicator_kmers'] = {}
            self.kmers[kmer]['non_indicator_kmers'] = {}
            self.kmers[kmer]['non_indicator_masks'] = {}
            n = 0
            for locus in self.kmers[kmer]['loci'].keys():
                l = len(list(filter(lambda x: x in self.kmers[kmer]['interest_kmers'], self.kmers[kmer]['loci'][locus]['kmers'])))
                m = len(list(filter(lambda x: x in self.kmers[kmer]['interest_masks'], self.kmers[kmer]['loci'][locus]['masks'])))
                if l != 0 or m != 0:
                    for ukmer in self.kmers[kmer]['loci'][locus]['kmers']:
                        self.kmers[kmer]['indicator_kmers'][ukmer] = True
                    n += 1
                else:
                    for ukmer in self.kmers[kmer]['loci'][locus]['kmers']:
                        self.kmers[kmer]['non_indicator_kmers'][ukmer] = True
                    self.kmers[kmer]['loci'].pop(locus, None)
            self.kmers[kmer].pop('interest_kmers', None)
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        self.kmers = {}
        i = 0
        for batch in self.load_output():
            print('adding batch', i)
            i += 1
            for kmer in batch:
                self.kmers[kmer] = batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        self.tracks = {}
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {}
                self.tracks[track][kmer] = self.kmers[kmer]
        print(len(self.tracks))
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountLociIndicatorKmersJob(map_reduce.FirstGenotypingJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'CountLociIndicatorKmersJob'
    _category = 'programming'
    _previous_job = FilterLociIndicatorKmersJob

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
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            self.kmers = json.load(json_file)
        with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        print(len(self.kmers), 'kmers')
        self.round_robin()

    def transform(self):
        c = config.Configuration()
        cpp_dir = os.path.join(os.path.dirname(__file__), '../../cpp')
        if c.accelerate:
            command = " " + str(self.index) + " " + self.get_current_job_directory() +  " " + c.fastq_file + " " + str(self.num_threads) + " " + "0"
            print(command)
            if c.debug:
                output = subprocess.call(os.path.join(cpp_dir, "counter_s.out") + command, shell = True)
            else:
                output = subprocess.call(os.path.join(cpp_dir, "counter.out") + command, shell = True)
            exit()
        else:
            with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'r') as json_file:
                self.kmers = json.load(json_file)
            for read, name, index in self.parse_fastq():
                self.process_read(read, name, index)
    
    def process_read(self, read, name, index):
        c = config.Configuration()
        slack = (c.read_length - c.ksize) / 2
        i = 0
        l = len(read)
        while i <= l - c.ksize:
            k = read[i : i + c.ksize]
            kmer = canonicalize(read[i: i + c.ksize])
            left = read[:i][-slack:]
            right = read[i + c.ksize:][:slack]
            if not kmer in self.kmers:
                i += 1
                continue
            found = False
            for locus in self.kmers[kmer]['loci']:
                for mask in self.kmers[kmer]['loci'][locus]['masks']:
                    if is_canonical_subsequence(mask, left) or is_canonical_subsequence(mask, right):
                        self.kmers[kmer]['count'] += 1
                        found = True
                        break
                if found:
                    debug_print(kmer)
                    debug_print(left)
                    debug_print(right)
                    debug_print(index)
                    debug_breakpoint()
                    break
            i += 1

    def compare(self):
        self.kmers = json.load(open(os.path.join(self.get_current_job_directory(), 'kmers.json')))
        self.py_kmers = {}
        self.tracks = {}
        print('loading python kmers')
        for kmer in self.kmers:
            for track in self.kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = True
                    print(track)
                    with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'r') as json_file:
                        kmers = json.load(json_file)
                        for kmer in kmers:
                            self.py_kmers[kmer] = kmers[kmer]['count']
        path = os.path.join(self.get_current_job_directory(), 'batch_0.json') 
        self.c_kmers = {}
        print(len(self.py_kmers))
        print('loading C kmers')
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                line = json_file.readline()
                while line:
                    kmer = line[:line.find(':')]
                    count = int(line[line.find(':') + 1:])
                    canon = canonicalize(kmer)
                    if not canon in self.c_kmers:
                        self.c_kmers[canon] = count
                    else:
                        self.c_kmers[canon] += count
                    line = json_file.readline()
        print(len(self.c_kmers))
        n = 0
        m = 0
        d = []
        s = 0
        print('comparing...')
        for kmer in self.py_kmers:
            if not kmer in self.c_kmers:
                n += 1
            else:
                if self.py_kmers[kmer] != self.c_kmers[kmer] / 2:
                    print(green(self.py_kmers[kmer]), blue(self.c_kmers[kmer] / 2))
                    a = abs(self.py_kmers[kmer] - self.c_kmers[kmer])
                    d.append(a)
                    s += a
                    m += 1
        visualizer.histogram(d, 'C-Python diff', self.get_current_job_directory(), 'diff', 'frequency') 
        print(n)
        print(m)
        print(float(s) / len(d))
        exit()

    def merge_accelerator_counts(self):
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
            n = 0
            with open (path, 'r') as json_file:
                line = json_file.readline()
                while line:
                    try:
                        kmer = line[:line.find(':')]
                        count = int(line[line.find(':') + 1:])
                        canon = canonicalize(kmer)
                        self.kmers[canon]['count'] += count / 2
                        line = json_file.readline()
                        n += 1
                    except Exception as e:
                        print(n, i, line)
                        print(e)
                        debug_breakpoint()
                        line = json_file.readline()
                        n += 1

    def reduce(self):
        c = config.Configuration()
        if c.accelerate:
            self.merge_accelerator_counts()
        else:
            self.kmers = self.merge_counts('total', 'doubt')
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
            with open(os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: 'indicator_kmers_' + track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class LociIndicatorKmersIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'LociIndicatorKmersIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = CountLociIndicatorKmersJob
    _kmer_type = 'non_unique_inner'

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = LociIndicatorKmersIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        c = config.Configuration()
        print(green(track_name))
        t = bed.track_from_name(track_name)
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            t = bed.track_from_name(track_name)
            kmers = json.load(json_file)
            if not kmers:
                print('no inner kmers found for', red(track_name))
                return None
            for kmer in kmers:
                self.lp_kmers[kmer] = {}
                self.lp_kmers[kmer]['type'] = self._kmer_type
                self.lp_kmers[kmer]['count'] = kmers[kmer]['count']# + k['doubt']
                self.lp_kmers[kmer]['doubt'] = kmers[kmer]['doubt']
                self.lp_kmers[kmer]['total'] = kmers[kmer]['total']
                self.lp_kmers[kmer]['reference'] = len(kmers[kmer]['loci'])
                self.lp_kmers[kmer]['reduction'] = kmers[kmer]['reference']
                self.lp_kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                self.lp_kmers[kmer]['loci'] = list(map(lambda l: l, kmers[kmer]['loci']))
        path = os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in kmers}, json_file, indent = 4, sort_keys = True)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenotypingConfidenceJob(map_reduce.BaseGenotypingJob):

    _name = 'GenotypingConfidenceJob'
    _category = 'programming'
    _previous_job = LociIndicatorKmersIntegerProgrammingJob
    _kmer_type = None

    @staticmethod
    def launch(**kwargs):
        job = GenotypingConfidenceJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        with open(os.path.join(self.get_previous_job_directory(), 'indexed_tracks.json'), 'r') as json_file:
            self.tracks = json.load(json_file)
        with open(os.path.join(self.get_previous_job_directory(), 'indexed_kmers.json'), 'r') as json_file:
            self.lp_kmers = json.load(json_file)
        self.round_robin(self.tracks)

    def transform(self, _, track):
        c = config.Configuration()
        print(cyan('probing', track))
        for r in [('00', 1.0), ('10', 0.5), ('11', 0.0)]:
            self.tracks[track]['objective_values_' + r[0]] = self.recalculate_lp_error(track, r[1])
            self.tracks[track]['objective_' + r[0]] = sum(self.tracks[track]['objective_values_' + r[0]])
        p = ['00', '10', '11']
        g = self.tracks[track]['lp_genotype']
        e = self.tracks[track]['objective_' + g]
        p.remove(g)
        m = min(p, key = lambda k: self.tracks[track]['objective_' + k])
        d = [a - b for a,b in zip(self.tracks[track]['objective_values_' + m], self.tracks[track]['objective_values_' + g])]
        n = len(d)
        s = statistics.std(d)
        t = sum(d) / (math.sqrt(n) * s) if s != 0 else float('Inf')
        t_value, p_value = scipy.stats.ttest_rel(self.tracks[track]['objective_values_' + m], self.tracks[track]['objective_values_' + g])
        self.tracks[track]['t_value'] = t_value
        self.tracks[track]['p_value'] = p_value
        #if math.isnan(t_value):
        #    print(yellow(track, p_value, t_value, t, d))
        #elif t_value == float('Inf'):
        #    print(magenta(track, p_value, t_value, t, d, self.tracks[track]['objective_values_' + m], self.tracks[track]['objective_values_' + g], c.coverage))
        #else:
        #    print(p_value, t_value, t)
        with open(os.path.join(self.get_current_job_directory(), 'confidence_' + track + '.json'), 'w') as json_file:
            json.dump(self.tracks[track], json_file, indent = 4)
        return None

    def recalculate_lp_error(self, track, c):
        errors = []
        for kmer in self.tracks[track]['kmers']:
            kmer = self.lp_kmers[kmer]
            r = 0
            for t in kmer['tracks']:
                if t == track:
                    r += kmer['coverage'] * kmer['tracks'][t] * (1.0 - 0.03) * c
                else:
                    r += kmer['coverage'] * kmer['tracks'][t] * (1.0 - 0.03) * self.tracks[t]['lp_rounding']
            errors.append(abs(kmer['count'] - kmer['coverage'] * kmer['residue'] - r))
        return errors

    def output_batch(self, batch):
        pass

    def reduce(self):
        for track in self.tracks:
            print('loading', track)
            with open(os.path.join(self.get_current_job_directory(), 'confidence_' + track + '.json'), 'r') as json_file:
                self.tracks[track] = json.load(json_file)
        self.export_confidence_scores()
        self.plot_confidence_scores()

    def export_confidence_scores(self):
        with open(os.path.join(self.get_current_job_directory(), 'confidence.bed'), 'w') as bed_file:
            for track in self.tracks:
                if str(self.tracks[track]['t_value']).find('inf') != -1:
                    self.tracks[track]['t_value'] = 100000
                if math.isnan(self.tracks[track]['t_value']):
                    continue
                t = bed.track_from_name(track)
                bed_file.write(t.chrom + '\t' +
                    str(t.begin) + '\t' +
                    str(t.end)   + '\t' +
                    str(self.tracks[track]['lp_value'])  + '\t' +
                    str(self.tracks[track]['lp_rounding'])  + '\t' +
                    str(self.tracks[track]['lp_genotype'])  + '\t' +
                    str(self.tracks[track]['actual_genotype'])  + '\t' +
                    str(self.tracks[track]['t_value'] if self.tracks[track]['p_value'] != float('Inf') else 1000) + '\t' +
                    str(self.tracks[track]['p_value'] if self.tracks[track]['t_value'] != float('Inf') else 1000) + '\t' +
                    str(len(self.tracks[track]['kmers'])) + '\n')

    def plot_confidence_scores(self):
        x = []
        p = []
        t = []
        for track in self.tracks:
            if math.isnan(self.tracks[track]['t_value']):
                continue
            x.append('Correct' if self.tracks[track]['lp_genotype'] == self.tracks[track]['actual_genotype'] else 'Wrong')
            p.append(self.tracks[track]['p_value'] if self.tracks[track]['p_value'] != float('Inf') else 1000)
            t.append(self.tracks[track]['t_value'] if self.tracks[track]['t_value'] != float('Inf') else 1000)
        visualizer.violin(x, p, 'p_value_distribution', self.get_current_job_directory(), 'Prediction', 'Confidence')
        visualizer.violin(x, t, 't_value_distribution', self.get_current_job_directory(), 'Prediction', 'Confidence')
        visualizer.histogram([t[i] for i in list(filter(lambda j: x[j] == 'Correct', range(0, len(x))))], 'confidence_score_correct_calls', self.get_current_job_directory(), 'Prediction', 'Confidence')
        visualizer.histogram([t[i] for i in list(filter(lambda j: x[j] == 'Wrong', range(0, len(x))))], 'confidence_score_wrong_calls', self.get_current_job_directory(), 'Prediction', 'Confidence')

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
