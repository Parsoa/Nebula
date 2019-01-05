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
    depth,
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
        #print(cyan(track))
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            t = bed.track_from_name(track_name)
            for kmer in kmers:
                #if kmer != 'CGGGTTCACGCCCTTCTCCCGCCTCAGCCTCC':
                #    continue
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

    def get_shared_masks(self, interest_kmers, kmers):
        l = []
        for x in kmers:
            for y in interest_kmers:
                if is_canonical_subsequence(x[4:-4], y):
                    l.append(x)
                    break
        return len(l)

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
                l = self.get_shared_masks(self.kmers[kmer]['interest_kmers'], self.kmers[kmer]['loci'][locus]['kmers'])
                #l = len(list(filter(lambda x: x in self.kmers[kmer]['interest_kmers'], self.kmers[kmer]['loci'][locus]['kmers'])))
                #m = len(list(filter(lambda x: x in self.kmers[kmer]['interest_masks'], self.kmers[kmer]['loci'][locus]['masks'])))
                #print(locus, l, list(filter(lambda x: x in self.kmers[kmer]['interest_kmers'], self.kmers[kmer]['loci'][locus]['kmers'])))
                if l != 0:
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
        visualizer.histogram(list(map(lambda kmer: len(self.kmers[kmer]['loci']), self.kmers)), 'number_of_loci', self.get_current_job_directory(), 'number of loci', 'number of kmers')
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
    _counter_mode = 0

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

    def reduce(self):
        c = config.Configuration()
        self.merge_accelerator_counts()
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

class AdjustCoverageGcContentJob(map_reduce.BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'AdjustCoverageGcContentJob'
    _category = 'programming'
    _previous_job = CountLociIndicatorKmersJob

    @staticmethod
    def launch(**kwargs):
        job = AdjustCoverageGcContentJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            self.kmers = json.load(json_file)
        print(len(self.kmers))
        gc_coverage_job = depth.ChromosomeGcContentEstimationJob()
        with open(os.path.join(gc_coverage_job.get_current_job_directory(), 'coverage.json'), 'r') as json_file:
            self.coverage = json.load(json_file)
        self.round_robin(self.kmers)
        self.batch_kmers = {}

    def transform(self, kmer, k):
        coverage = []
        for locus in kmer['loci']:
            gc = calculate_gc_content(kmer['loci'][locus]['seq']['all'])
            coverage.append(self.coverage[gc])
        kmer['coverage'] = statistics.mean(coverage)
        self.batch_kmers[k] = kmer
        return None

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.batch_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def reduce(self):
        self.kmers = {}
        for batch in self.load_output():
            for kmer in batch:
                self.kmers[kmer] = batch[kmer]
        print(len(self.kmers))
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
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
    _previous_job = AdjustCoverageGcContentJob
    _kmer_type = 'inner'

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
            if len(kmers) == 0:
                print('no inner kmers found for', red(track_name))
                return None
            #kmers = {kmer: kmer[kmer] for kmer in list(filter(lambda kmer: len(kmers[kmer]['tracks']) == 1, kmers))}
            unique_kmers = list(filter(lambda kmer: len(kmers[kmer]['loci']) == 1, kmers))
            if len(unique_kmers) > 3:
                kmers = {kmer: kmers[kmer] for kmer in unique_kmers}
            for kmer in kmers:
                self.lp_kmers[kmer] = {}
                self.lp_kmers[kmer]['type'] = self._kmer_type
                self.lp_kmers[kmer]['count'] = kmers[kmer]['count']# if kmers[kmer]['reference'] == 1 else kmers[kmer]['total']
                self.lp_kmers[kmer]['doubt'] = kmers[kmer]['doubt']
                self.lp_kmers[kmer]['total'] = kmers[kmer]['total']
                self.lp_kmers[kmer]['coverage'] = kmers[kmer]['coverage']
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
# ============================================================================================================================ #
# ============================================================================================================================ #

class LociIndicatorKmersFractionalProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'LociIndicatorKmersFractionalProgrammingJob'
    _category = 'programming'
    _previous_job = AdjustCoverageGcContentJob
    _kmer_type = 'inner'

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = LociIndicatorKmersFractionalProgrammingJob(**kwargs)
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
            if len(kmers) == 0:
                print('no inner kmers found for', red(track_name))
                return None
            _kmers = {kmer: kmers[kmer] for kmer in list(filter(lambda kmer: len(kmers[kmer]['tracks']) == 1, kmers))}
            unique_kmers = list(filter(lambda kmer: len(_kmers[kmer]['loci']) == 1, _kmers))
            non_unique_kmers = list(filter(lambda kmer: len(_kmers[kmer]['loci']) != 1, _kmers))
            # let's limit the number of kmers to 50 to help with integer programming runtime
            if len(unique_kmers) >= 5:
                kmers = {kmer: _kmers[kmer] for kmer in unique_kmers[:min(50, len(unique_kmers))]}
            else:
                kmers = {kmer: _kmers[kmer] for kmer in unique_kmers}
                for kmer in non_unique_kmers:
                    kmers[kmer] = _kmers[kmer]
                    if len(kmers) > 50:
                        break
            for kmer in kmers:
                self.lp_kmers[kmer] = {}
                self.lp_kmers[kmer]['type'] = self._kmer_type
                self.lp_kmers[kmer]['count'] = kmers[kmer]['count']# if kmers[kmer]['reference'] == 1 else kmers[kmer]['total']
                self.lp_kmers[kmer]['doubt'] = kmers[kmer]['doubt']
                self.lp_kmers[kmer]['total'] = kmers[kmer]['total']
                self.lp_kmers[kmer]['coverage'] = kmers[kmer]['coverage']
                self.lp_kmers[kmer]['reference'] = len(kmers[kmer]['loci'])
                self.lp_kmers[kmer]['reduction'] = kmers[kmer]['reference']
                self.lp_kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                self.lp_kmers[kmer]['loci'] = list(map(lambda l: l, kmers[kmer]['loci']))
        path = os.path.join(self.get_current_job_directory(), 'indicator_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in kmers}, json_file, indent = 4, sort_keys = True)
        return path

    def solve(self):
        helper = LociIndicatorKmersFractionalProgrammingJob.FractionalProgramHelper()
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
                g = self.round_genotype(c / w)
                bed_file.write(t.chrom + '\t' +
                    str(t.begin) + '\t' +
                    str(t.end) + '\t' +
                    str(g[1]) + '\t' +
                    str(1.0 / w) + '\t' +
                    str(c / w) + '\t' +
                    str(len(solution) / 2) + '\n')
                with open(os.path.join(self.get_current_job_directory(), 'solution_' + track + '.json'), 'w') as json_file:
                    json.dump({'variables': [s / w if s != w else 1.0 / s for s in solution]}, json_file, indent = 4, sort_keys = True)

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    class FractionalProgramHelper(programming.IntegerProgrammingJob):

        _name = 'LociIndicatorKmersFractionalProgrammingJob'
        _category = 'programming'
        _previous_job = AdjustCoverageGcContentJob
        _kmer_type = 'inner'

        def load_inputs(self):
            self.round_robin(self.tracks)

        def transform(self, _, track):
            print('solving for', cyan(track))
            problem = self.generate_linear_program(track)
            print('program_' + track + '.lp')
            problem.write(os.path.join(self.get_current_job_directory(), 'program_' + str(track) + '.lp'))
            problem.solve()
            solution = problem.solution.get_values()
            with open(os.path.join(self.get_current_job_directory(), 'solution_' + track + '.json'), 'w') as json_file:
                json.dump({'variables': solution}, json_file, indent = 4, sort_keys = True)
            return None

        def reduce(self):
            pass

        def output_batch(self, batch):
            exit()

        def generate_linear_program(self, track):
            print('generating linear program')
            c = config.Configuration()
            globals()['cplex'] = __import__('cplex')
            problem = cplex.Cplex()
            problem.objective.set_sense(problem.objective.sense.minimize)
            lp_kmers = list(filter(lambda kmer: track in kmer['tracks'], self.lp_kmers))
            # the coverage of each event
            tokens = track.split('_')
            problem.variables.add(names = ['y_c' + tokens[1]])
            # the real-valued error parameter for inner_kmer
            problem.variables.add(names = ['y_e' + str(index) for index, kmer in enumerate(lp_kmers)],
                lb = [-100000 for index, kmer in enumerate(lp_kmers)],
            )
            # absolute value of the inner_kmer error parameter
            problem.variables.add(names = ['y_l' + str(index) for index, kmer in enumerate(lp_kmers)],
                obj = [1.0] * len(lp_kmers),
            )
            # indicator variables
            problem.variables.add(names = ['y_i' + str(index) for index, kmer in enumerate(lp_kmers)], ub = [1.0] * len(lp_kmers), types = [problem.variables.type.integer] * len(lp_kmers))
            # t variable
            problem.variables.add(names = ['t'], lb = [-100000])
            # range constraints
            ind = [0, 1 + 3 * len(lp_kmers)]
            val = [1, -1]
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [0],
                senses = ['L']
            )
            # constraints
            ind = [0, 1 + 3 * len(lp_kmers)]
            val = [17, 3]
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [1],
                senses = ['E']
            )
            # indicators
            ind = [1 + 2 * len(lp_kmers) + index for index, kmer in enumerate(lp_kmers)]
            val = [1.0] * len(lp_kmers)
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [0.1 * len(lp_kmers)],
                senses = ['L']
            )
            n = 0
            start = time.time()
            for index, kmer in enumerate(lp_kmers):
                ind = [0] 
                ind.append(1 + index)
                ind.append(1 + 3 * len(lp_kmers))
                val = list(map(lambda track: kmer['coefficient'] * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                val.append(1.0)
                val.append(-1 * (kmer['count'] - kmer['coverage'] * kmer['residue']))
                #problem.linear_constraints.add(
                #    lin_expr = [cplex.SparsePair(
                #        ind = ind,
                #        val = val,
                #    )],
                #    rhs = [0],
                #    senses = ['E']
                #)
                problem.indicator_constraints.add(
                    indvar = 1 + 2 * len(lp_kmers) + index,
                    complemented = 1,
                    rhs = 0,
                    lin_expr = cplex.SparsePair(
                        ind = ind,
                        val = val,
                    ),
                    sense = 'E'
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

# ============================================================================================================================ #
# ============================================================================================================================ #
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
        with open(os.path.join(self.get_previous_job_directory(), 'tracks.json'), 'r') as json_file:
            self.tracks = json.load(json_file)
        with open(os.path.join(self.get_previous_job_directory(), 'lp_kmers.json'), 'r') as json_file:
            self.lp_kmers = json.load(json_file)
        with open(os.path.join(self.get_previous_job_directory(), 'solution.json'), 'r') as json_file:
            self.solution = json.load(json_file)['variables']
        #self.errors = sorted(self.solution[len(self.tracks):])
        #self.correct_errors = sorted(list(reduce(lambda x, y : x + y, [] + list(map(lambda track: [self.errors[i] for i in self.tracks[track]['kmers']], list(filter(lambda t: self.tracks[t]['lp_genotype'] == self.tracks[t]['actual_genotype'], self.tracks)))))))
        #self.wrong_errors = sorted(list(reduce(lambda x, y: x + y, [] + list(map(lambda track: [self.errors[i] for i in self.tracks[track]['kmers']], list(filter(lambda t: self.tracks[t]['lp_genotype'] != self.tracks[t]['actual_genotype'], self.tracks)))))))
        #self.estimate_distribution(self.correct_errors)
        #self.estimate_distribution(self.wrong_errors)
        #self.estimate_distribution(self.errors)
        #visualizer.histogram(self.correct_errors, 'correct_errors', self.get_current_job_directory(), 'error', 'frequency')
        #visualizer.histogram(self.wrong_errors, 'wrong_errors', self.get_current_job_directory(), 'error', 'frequency')
        #visualizer.histogram(self.errors, 'errors', self.get_current_job_directory(), 'error', 'frequency')
        #exit()
        self.round_robin(self.tracks)

    def estimate_distribution(self, x):
        std = statistics.std(x)
        mean = statistics.mean(x)
        median = statistics.median(x)
        print(mean, std, median)
        x = list(filter(lambda x: x < 3 * mean, x))
        std = statistics.std(x)
        mean = statistics.mean(x)
        median = statistics.median(x)
        print(mean, std, median)
        return
        x = list(filter(lambda x: x < 3 * mean, x))
        std = statistics.std(x)
        mean = statistics.mean(x)
        print(mean, std)

    def transform(self, _, track):
        print(cyan(track))
        l = len(self.lp_kmers)
        current = 0
        n = len(self.tracks[track]['kmers'])
        self.tracks[track]['bootstrap'] = {'00': [], '10': [], '11': [], 'genotypes': []}
        kmers = list(filter(lambda k: len(self.lp_kmers[k]['tracks']) == 1, self.tracks[track]['kmers']))
        kmers = [self.lp_kmers[i] for i in self.tracks[track]['kmers']]
        for i in range(n):
            kmers[i]['coefficient'] = 0
            if len(kmers[i]['tracks']) > 1:
                for t in kmers[i]['tracks']:
                    if t != track:
                        kmers[i]['count'] -= self.tracks[t]['lp_value'] * 0.97 * kmers[i]['coverage'] * kmers[i]['tracks'][t]
        for i in range(n):
            print(green('iterarion', i))
            kmers[i]['coefficient'] = 1
            problem = self.generate_linear_program(kmers, track)
            problem.solve()
            self.solution = problem.solution.get_values()
            s = int(round(2 * self.solution[0]))
            g = '00' if s == 2 else '10' if s == 1 else '11'
            self.tracks[track]['bootstrap']['genotypes'].append(g)
            self.tracks[track]['bootstrap'][g].append(i)
        with open(os.path.join(self.get_current_job_directory(), 'bootstrap_' + track + '.json'), 'w') as json_file:
            json.dump(self.tracks[track], json_file, indent = 4)
        return None

    def generate_linear_program(self, kmers, track):
        c = config.Configuration()
        globals()['cplex'] = __import__('cplex')
        problem = cplex.Cplex()
        problem.set_log_stream(None)
        problem.set_error_stream(None)
        problem.set_warning_stream(None)
        problem.set_results_stream(None)
        problem.objective.set_sense(problem.objective.sense.minimize)
        # the coverage of each event
        problem.variables.add(names = ['t'],
            ub = [1.0],
        )
        # the real-valued error parameter for inner_kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(kmers)],
            ub = [kmer['coefficient'] * (kmer['count'] - kmer['coverage'] * kmer['residue']) for kmer in kmers],
            lb = [kmer['coefficient'] * (kmer['count'] - kmer['coverage'] * kmer['residue'] - kmer['coverage'] * sum(kmer['tracks'][track] for track in kmer['tracks'])) for kmer in kmers],
        )
        # absolute value of the inner_kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(kmers)],
            obj = [1.0] * len(kmers),
        )
        #self.add_snp_linear_constraint(problem)
        n = 0
        start = time.time()
        offset = len(self.tracks) + 2 * len(kmers)
        for index, kmer in enumerate(kmers):
            if kmer['coefficient'] == 0:
                continue
            # TxR + E = C - 
            #ref = kmer['reference']
            ind = [0]
            #ind += [offset + i for i in range(0, ref)] #SNP
            ind.append(1 + index) # Objective
            val = [kmer['coefficient'] * kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03)] #Coverage corrected for errors
            #val += [-kmer['coverage']] * ref #SNP
            val.append(1.0) #Objective
            #offset += ref
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = val,
                )],
                rhs = [kmer['coefficient'] * (kmer['count'] - kmer['coverage'] * kmer['residue'])],
                senses = ['E']
            )
            self.add_error_absolute_value_constraints(problem, index, kmers)
            n = n + 1
            if n % 1000 == 0:
                t = time.time()
                p = float(n) / len(kmers)
                eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
                #print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
        return problem

    def add_error_absolute_value_constraints(self, problem, index, kmers):
        problem.linear_constraints.add(
            lin_expr = [cplex.SparsePair(
                ind = [1 + len(kmers) + index, 1 + index],
                val = [1.0, 1.0],
            )],
            rhs = [0],
            senses = ['G']
        )
        problem.linear_constraints.add(
            lin_expr = [cplex.SparsePair(
                ind = [1 + len(kmers) + index, 1 + index],
                val = [1.0, -1.0],
            )],
            rhs = [0],
            senses = ['G']
        )

    def reduce(self):
        remove = []
        for track in self.tracks:
            path = os.path.join(self.get_current_job_directory(), 'bootstrap_' + track + '.json')
            try:
                with open(path, 'r') as json_file:
                    self.tracks[track] = json.load(json_file)
            except:
                remove.append(track)
        for track in remove:
            self.tracks.pop(track, None)
        with open(os.path.join(self.get_current_job_directory(), 'bootstrap.bed'), 'w') as bed_file:
            for track in self.tracks:
                t = bed.track_from_name(track)
                bed_file.write(t.chrom + '\t' +
                    str(t.begin) + '\t' +
                    str(t.end)   + '\t' +
                    str(self.tracks[track]['lp_value'])  + '\t' +
                    str(self.tracks[track]['lp_rounding'])  + '\t' +
                    str(self.tracks[track]['lp_genotype'])  + '\t' +
                    str(self.tracks[track]['actual_genotype']) + '\t' +
                    str(len(self.tracks[track]['kmers'])) + '\t' +
                    str(self.tracks[track]['confidence_score']) + '\t' +
                    compactify(self.tracks[track]['bootstrap']['00'])  + '\t' +
                    compactify(self.tracks[track]['bootstrap']['10'])  + '\t' +
                    compactify(self.tracks[track]['bootstrap']['11'])  + '\t' +
                    compactify(self.tracks[track]['bootstrap']['genotypes'])  + '\n')
        self.plot_rounding_gap()

    def plot_rounding_gap(self):
        x = []
        s = []
        for track in self.tracks:
            x.append('Correct' if self.tracks[track]['lp_genotype'] == self.tracks[track]['actual_genotype'] else 'Wrong')
            if len(self.tracks[track]['bootstrap']['genotypes']) == 1:
                s.append(-1)
            else:
                g = self.tracks[track]['bootstrap']['genotypes']
                s.append(sum(map(lambda i: 1 if g[i] != g[i - 1] else 0, range(1, len(g))))) 
        visualizer.violin(x, s, 'rounding_gap_distribution', self.get_current_job_directory(), 'Prediction', 'Confidence')

    #def transform(self, _, track):
    #    c = config.Configuration()
    #    std = 5
    #    print(cyan('probing', track))
    #    stats = {'00': {}, '10': {}, '11': {}}
    #    for r in [('00', 1.0), ('10', 0.5), ('11', 0.0)]:
    #        self.tracks[track]['objective_values_' + r[0]] = self.recalculate_lp_error(track, r[1])
    #        self.tracks[track]['objective_' + r[0]] = sum(self.tracks[track]['objective_values_' + r[0]])
    #        mean = statistics.mean(self.tracks[track]['objective_values_' + r[0]])
    #        N = statistics.NormalDistribution(mean, std)
    #        self.tracks[track]['likelihood_' + r[0]] = reduce(lambda x, y: x * N.pmf(int(y)), [1] + self.tracks[track]['objective_values_' + r[0]]) 
    #    p = ['00', '10', '11']
    #    g = self.tracks[track]['lp_genotype']
    #    e = self.tracks[track]['objective_' + g]
    #    p.remove(g)
    #    m = max(p, key = lambda k: self.tracks[track]['likelihood_' + k])
    #    N = statistics.NormalDistribution(0, std)
    #    self.tracks[track]['confidence_score'] = reduce(lambda x, y: x * N.pmf(int(y)), [1] + self.tracks[track]['objective_values_' + g]) / self.tracks[track]['likelihood_' + m]
    #    with open(os.path.join(self.get_current_job_directory(), 'confidence_' + track + '.json'), 'w') as json_file:
    #        json.dump(self.tracks[track], json_file, indent = 4)
    #    return None

    #def recalculate_lp_error(self, track):
    #    errors = []
    #    for kmer in self.tracks[track]['kmers']:
    #        kmer = self.lp_kmers[kmer]
    #        r = 0
    #        for t in kmer['tracks']:
    #            if t == track:
    #                r += kmer['coverage'] * kmer['tracks'][t] * (1.0 - 0.03) * self.tracks[t]['lp_rounding']
    #            else:
    #                r += kmer['coverage'] * kmer['tracks'][t] * (1.0 - 0.03) * self.tracks[t]['lp_rounding']
    #        errors.append(abs(kmer['count'] - kmer['coverage'] * kmer['residue'] - r))
    #    return errors

    #def output_batch(self, batch):
    #    pass

    #def reduce(self):
    #    for track in self.tracks:
    #        print('loading', track)
    #        with open(os.path.join(self.get_current_job_directory(), 'confidence_' + track + '.json'), 'r') as json_file:
    #            self.tracks[track] = json.load(json_file)
    #    self.export_confidence_scores()
    #    self.plot_confidence_scores()

    #def export_confidence_scores(self):
    #    with open(os.path.join(self.get_current_job_directory(), 'confidence.bed'), 'w') as bed_file:
    #        for track in self.tracks:
    #            t = bed.track_from_name(track)
    #            bed_file.write(t.chrom + '\t' +
    #                str(t.begin) + '\t' +
    #                str(t.end)   + '\t' +
    #                str(self.tracks[track]['lp_value'])  + '\t' +
    #                str(self.tracks[track]['lp_rounding'])  + '\t' +
    #                str(self.tracks[track]['lp_genotype'])  + '\t' +
    #                str(self.tracks[track]['actual_genotype']) + '\t' +
    #                str(self.tracks[track]['standard_error'])  + '\t' +
    #                str(len(self.tracks[track]['kmers'])) + '\n')

    #def plot_confidence_scores(self):
    #    x = []
    #    e = []
    #    for track in self.tracks:
    #        x.append('Correct' if self.tracks[track]['lp_genotype'] == self.tracks[track]['actual_genotype'] else 'Wrong')
    #        e.append(self.tracks[track]['standard_error'])
    #    visualizer.violin(x, e, 'standard_error_distribution', self.get_current_job_directory(), 'Prediction', 'Confidence')

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
