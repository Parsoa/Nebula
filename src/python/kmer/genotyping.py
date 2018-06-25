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
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import plotly.offline as plotly
import plotly.graph_objs as graph_objs

import numpy

# ============================================================================================================================ #
# MapReduce job to genotype a genome
# We need to have generated a counttable for the genome we are interested in beforehand
# How to get this far?
# Step 1: get the break points all of which's kmers are found in CHM1, these are the possible candidates for the structural
# variation's boundaries -> BreakPointJob
# Step 2: Find a set of novel kmers for each break point that can be used to indentify it. khmer never underestimates counts so
# if a kmer comes with a count of zero in reference genome, we can be sure that it is really novel -> NovelKmerJob
# Step 3: khmer may report oversetimated counts for these break points so we need to count them exactly again. This is necessary
# for a reliable likelihood model -> CountKmersExactJob
# Step 4: With exact kmer counts available, we can find the most likely break points for each event in our library -> MostLikelyBreakPointsJob
# Step 6: For every given genome we will need to calculate the coverage depth ->
# Step 5: Given a sample genome, try to genotype the structural variations using the likelihood model and signatures gathered
# above -> GenotypingJob (this one)
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenotypingJob(map_reduce.BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = GenotypingJob(job_name = 'Genotyping_', previous_job_name = 'MostLikelyBreakPointsJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare(self):
        c = config.Configuration()
        self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])

    def transform(self, track, track_name):
        c = config.Configuration()
        likelihood = {}
        distribution = {
            '(1, 1)': statistics.NormalDistribution(mean = c.coverage, std = c.std),
            '(1, 0)': statistics.NormalDistribution(mean = c.coverage / 2, std = c.std),
            '(0, 0)': statistics.NormalDistribution(mean = 0, std = c.std)
        }
        inner_distribution = {
            '(0, 0)': statistics.NormalDistribution(mean = c.coverage, std = c.std),
            '(1, 0)': statistics.NormalDistribution(mean = c.coverage / 2, std = c.std),
            '(1, 1)': statistics.NormalDistribution(mean = 0, std = c.std)
        }
        # all kmers remaining in the track are to be considered part of the signature, no matter the number of breakpoints
        with open(track, 'r') as track_file:
            break_points = json.load(track_file)
        for break_point in break_points:
            if len(break_points[break_point]['left_local_unique_kmers']) != 0:
                return None
            if len(break_points[break_point]['right_local_unique_kmers']) != 0:
                return None
            likelihood[break_point] = {}
            likelihood[break_point]['novel_kmers'] = {}
            likelihood[break_point]['inner_kmers'] = {}
            likelihood[break_point]['2std'] = 0
            likelihood[break_point]['likelihood'] = {}
            likelihood[break_point]['likelihood']['novel'] = {}
            likelihood[break_point]['likelihood']['inner'] = {}
            for zyg in distribution:
                likelihood[break_point]['likelihood']['inner'][zyg] = 0
                likelihood[break_point]['likelihood']['novel'][zyg] = 0
            for kmer in break_points[break_point]['novel_kmers']:
                count = self.counts_provider.get_kmer_count(kmer)
                if count > c.coverage + 3 * c.std:
                    continue
                likelihood[break_point]['novel_kmers'][kmer] = count
                for zyg in distribution:
                    likelihood[break_point]['likelihood']['novel'][zyg] += distribution[zyg].log_pmf(count)
            for kmer in break_points[break_point]['inner_kmers']:
                count = self.counts_provider.get_kmer_count(kmer)
                if count > c.coverage + 3 * c.std:
                    continue
                likelihood[break_point]['inner_kmers'][kmer] = count
                for zyg in distribution:
                    likelihood[break_point]['likelihood']['inner'][zyg] += inner_distribution[zyg].log_pmf(count)
            for zyg in distribution:
                likelihood[break_point]['likelihood']['inner'][zyg] -= len(likelihood[break_point]['inner_kmers']) * distribution[zyg].log_pmf(distribution[zyg].mean)
                likelihood[break_point]['likelihood']['novel'][zyg] -= len(likelihood[break_point]['novel_kmers']) * distribution[zyg].log_pmf(distribution[zyg].mean)
                likelihood[break_point]['likelihood']['inner'][zyg] = abs(likelihood[break_point]['likelihood']['inner'][zyg])
                likelihood[break_point]['likelihood']['novel'][zyg] = abs(likelihood[break_point]['likelihood']['novel'][zyg])
            inner_choice, inner_cost = min(likelihood[break_point]['likelihood']['inner'].items(), key = operator.itemgetter(1))
            novel_choice, novel_cost = min(likelihood[break_point]['likelihood']['novel'].items(), key = operator.itemgetter(1))
            likelihood[break_point]['genotype'] = {}
            # too few kmers in to do anything useful
            #if len(likelihood[break_point]['inner_kmers']) <= 5 and len(likelihood[break_point]['novel_kmers']) <= 5:
                #choice = '(2, 2)'
            # enough novel kmers but not enough inner kmers
            if len(likelihood[break_point]['inner_kmers']) <= 5:
                inner_choice = '(3, 3)'
                choice = inner_choice
            # enough inner kmers but not enough novel kmers
            elif len(likelihood[break_point]['novel_kmers']) <= 5:
                novel_choice = '(4, 4)'
                choice = inner_choice
            else:
                inner_cost = inner_cost / len(likelihood[break_point]['inner_kmers'])
                novel_cost = novel_cost / len(likelihood[break_point]['novel_kmers'])
                choice = inner_choice# if inner_cost < novel_cost else novel_choice
            likelihood[break_point]['genotype']['inner'] = inner_choice
            likelihood[break_point]['genotype']['novel'] = novel_choice
            likelihood[break_point]['genotype']['consensus'] = choice
        path = os.path.join(self.get_current_job_directory(), 'genotype_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(likelihood, json_file, sort_keys = True, indent = 4)
        return path

    def reduce(self):
        c = config.Configuration()
        output = map_reduce.Job.reduce(self)
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            #bed_file.write('chrom\tstart\tend\tkmers\tgenotype\t0,0\t1,0\t1,1\tcorrect\n')
            for track in output:
                with open(output[track], 'r') as track_file:
                    break_points = json.load(track_file)
                    chrom = track.split('_')[0]
                    begin = track.split('_')[1]
                    end = track.split('_')[2]
                    for break_point in break_points:
                        inner_kmers = len(break_points[break_point]['inner_kmers'])
                        novel_kmers = len(break_points[break_point]['novel_kmers'])
                        inner_zyg00 = break_points[break_point]['likelihood']['inner']['(0, 0)']
                        inner_zyg10 = break_points[break_point]['likelihood']['inner']['(1, 0)']
                        inner_zyg11 = break_points[break_point]['likelihood']['inner']['(1, 1)']
                        novel_zyg00 = break_points[break_point]['likelihood']['novel']['(0, 0)']
                        novel_zyg10 = break_points[break_point]['likelihood']['novel']['(1, 0)']
                        novel_zyg11 = break_points[break_point]['likelihood']['novel']['(1, 1)']
                        inner = break_points[break_point]['genotype']['inner']
                        novel = break_points[break_point]['genotype']['novel']
                        consensus = break_points[break_point]['genotype']['consensus']
                        bed_file.write(chrom + '\t' + begin + '\t' + end + '\t' + str(novel_kmers) + '\t' + str(inner_kmers) + '\t' +
                                str(inner_zyg00 + novel_zyg00) + '\t' +
                                str(inner_zyg10 + novel_zyg10) + '\t' +
                                str(inner_zyg11 + novel_zyg11) + '\t' +
                                str(consensus) + '\t' +
                                str(inner) + '\t' +
                                str(novel) + '\t' +
                                '\n')
        return output

    def plot(self, tracks):
        inner = []
        novel = []
        for track in tracks:
            with open(tracks[track], 'r') as track_file:
                break_points = json.load(track_file)
                for break_point in break_points:
                    inner.append(len(break_points[break_point]['inner_kmers']))
                    novel.append(len(break_points[break_point]['novel_kmers']))
        visualizer.histogram(inner, 'num_unique_inner_kmers', self.get_current_job_directory(), x_label = 'number of unique inner kmers', y_label = 'number of events')
        visualizer.histogram(inner, 'num_novel_kmers', self.get_current_job_directory(), x_label = 'number of novel kmers', y_label = 'number of events')
        #
        self.plot_novel_kmer_count_in_sample(tracks)

    def plot_novel_kmer_count_in_sample(self, tracks):
        novel = []
        for track in tracks:
            with open(tracks[track], 'r') as track_file:
                break_points = json.load(track_file)
                for break_point in break_points:
                    for kmer in break_points[break_point]['novel_kmers']:
                        novel.append(break_points[break_point]['novel_kmers'][kmer])
        visualizer.histogram(novel, 'novel_kmers_in_sample', self.get_current_job_directory(), x_label = 'number of times novel kmers appears in sample', y_label = 'number of kmers')

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        d = 'simulation' if c.simulation else 'output'
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../' + d + '/' + bed_file_name + '/' + str(c.ksize) + '/', self.previous_job_name[:-1]))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class LocalUniqueKmersCountingJob(map_reduce.BaseGenotypingJob):

    def get_next_read(self):
        name = self.fastq_file.readline().strip()
        read = self.fastq_file.readline().strip()
        self.fastq_file.readline().strip()
        self.fastq_file.readline().strip()
        return (read, name)

    def parse_paired_end_fastq(self):
        line = self.fastq_file.readline().strip()
        if len(line) == 0:
            line = self.fastq_file.readline().strip()
        while True:
            if line[0] == '@' and line.endswith('/1'):
                right = (self.fastq_file.readline().strip(), line)
                break
            line = self.fastq_file.readline().strip()
        self.fastq_file.readline().strip()
        self.fastq_file.readline().strip()
        n = 0
        m = 0
        t = time.time()
        LEFT = 1
        RIGHT = 0
        state = LEFT
        while True:
            read = self.get_next_read()
            if read[0] == '':
                break
            if state == RIGHT:
                right = read
                n += 4
                if right[1].endswith('/1'):
                    state = LEFT
            else:
                left = read
                n += 4
                if left[1].endswith('/2'):
                    if right[1][0: -2] == left[1][0: -2]:
                        #print('read', left[1], right[1])
                        yield right[0], left[0]
                        state = RIGHT
                        if self.fastq_file.tell() >= (self.index + 1) * self.fastq_file_chunk_size:
                            break
                else:
                    right = left
                    state = LEFT
            if n >= 100000:
                n = 0
                m += 1
                c = self.fastq_file.tell() - self.index * self.fastq_file_chunk_size
                s = time.time()
                p = c / float(self.fastq_file_chunk_size)
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:2d}'.format(self.index), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = LocalUniqueKmersCountingJob(job_name = 'Genotyping_', previous_job_name = 'MostLikelyBreakPointsJob_', batch_file_prefix = 'experimental_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    # needs to know the most likely breakpoint and its kmers, MostLikelyBreakPointsJob includes that information
    def load_inputs(self):
        # each event is a structural variation with its most likely breakpoints
        self.left_local_unique_kmers = {}
        self.right_local_unique_kmers = {}
        self.tracks = {}
        tracks = self.load_previous_job_results()
        for track in tracks:
            break_points = json.load(open(tracks[track], 'r'))
            for break_point in break_points:
                # only run for events that don't have standalone inner kmers
                if len(break_points[break_point]['left_local_unique_kmers']) != 0:
                    self.tracks[track] = break_points
        duplicates = {}
        for track in self.tracks:
            for break_point in self.tracks[track]:
                self.tracks[track][break_point].pop('likelihood', None)
                self.tracks[track][break_point]['inner_kmers'] = {kmer: 0 for kmer in self.tracks[track][break_point]['inner_kmers']}
                for kmer in self.tracks[track][break_point]['left_local_unique_kmers']:
                    if kmer in self.left_local_unique_kmers:
                        print(red('duplicate'))
                        duplicates[kmer] = 0
                    self.left_local_unique_kmers[kmer] = {
                        'count': 0,
                        'track': track,
                        'break_point': break_point,
                        'original': kmer
                    }
                    self.left_local_unique_kmers[reverse_complement(kmer)] = {
                        'count': 0,
                        'track': track,
                        'break_point': break_point,
                        'original': kmer
                    }
                for kmer in self.tracks[track][break_point]['right_local_unique_kmers']:
                    if kmer in self.right_local_unique_kmers:
                        print(red('duplicate'))
                        duplicates[kmer] = 0
                    self.right_local_unique_kmers[kmer] = {
                        'count': 0,
                        'track': track,
                        'break_point': break_point,
                        'original': kmer
                    }
                    self.right_local_unique_kmers[reverse_complement(kmer)] = {
                        'count': 0,
                        'track': track,
                        'break_point': break_point,
                        'original': kmer
                    }
            self.num_threads = c.max_threads
            for i in range(self.num_threads):
                self.batch[i] = {}
        print('duplicates:', green(len(duplicates)), 'total:', blue(len(self.left_local_unique_kmers) + len(self.right_local_unique_kmers)))

    def run_batch(self, batch):
        c = config.Configuration()
        self.fastq_file = open(c.fastq_file, 'r')
        self.fastq_file_chunk_size = math.ceil(os.path.getsize(self.fastq_file.name) / float(self.num_threads))
        self.fastq_file.seek(self.index * self.fastq_file_chunk_size, 0)
        # this forked process will exit at the end of the following function call
        self.transform()
        self.output_batch(self.tracks)
    
    def transform(self):
        for left, right in self.parse_paired_end_fastq():
            #print(blue(left), green(right))
            self.update_counts(right = right, left = left, right_kmers = self.right_local_unique_kmers, left_kmers = self.left_local_unique_kmers)
            self.update_counts(right = right, left = left, right_kmers = self.left_local_unique_kmers, left_kmers = self.right_local_unique_kmers)
            self.update_counts(right = left, left = right, right_kmers = self.right_local_unique_kmers, left_kmers = self.left_local_unique_kmers)
            self.update_counts(right = left, left = right, right_kmers = self.left_local_unique_kmers, left_kmers = self.right_local_unique_kmers)

    def update_counts(self, right, left, right_kmers, left_kmers):
        found = {}
        for kmer in extract_kmers(c.ksize, left):
            #kmer = find_kmer(k, left_kmers)
            if kmer in left_kmers:
                track = left_kmers[kmer]['track']
                break_point = left_kmers[kmer]['break_point']
                original_kmer = left_kmers[kmer]['original']
                #if original_kmer in self.tracks[track][break_point]['left_local_unique_kmers']:
                #    self.tracks[track][break_point]['left_local_unique_kmers'][original_kmer] += 1
                #else:
                #    self.tracks[track][break_point]['right_local_unique_kmers'][original_kmer] += 1
                for right_kmer in extract_kmers(c.ksize, right):
                    if right_kmer in right_kmers:
                        if right_kmers[right_kmer]['track'] == left_kmers[kmer]['track']:
                            if original_kmer in self.tracks[track][break_point]['left_local_unique_kmers']:
                                self.tracks[track][break_point]['left_local_unique_kmers'][original_kmer] += 1
                            else:
                                self.tracks[track][break_point]['right_local_unique_kmers'][original_kmer] += 1
                            break
                if not track in found:
                    for inner_kmer in extract_kmers(c.ksize, right):
                        if inner_kmer in self.tracks[track][break_point]['inner_kmers']:
                            found[track] = True
                            self.tracks[track][break_point]['inner_kmers'][inner_kmer] += 1

    def reduce(self):
        c = config.Configuration()
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + '0.json')
        with open(path, 'r') as json_file:
            output = json.load(json_file)
        for i in range(1, self.num_threads):
            tracks = self.load_output_batch(i)
            for track in tracks:
                for break_point in tracks[track]:
                    for kmer in tracks[track][break_point]['inner_kmers']:
                        output[track][break_point]['inner_kmers'][kmer] += tracks[track][break_point]['inner_kmers'][kmer]
                    for kmer in tracks[track][break_point]['left_local_unique_kmers']:
                        output[track][break_point]['left_local_unique_kmers'][kmer] += tracks[track][break_point]['left_local_unique_kmers'][kmer]
                    for kmer in tracks[track][break_point]['right_local_unique_kmers']:
                        output[track][break_point]['right_local_unique_kmers'][kmer] += tracks[track][break_point]['right_local_unique_kmers'][kmer]
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)

    def count_inner_kmers(self, read):
        kmers = extract_kmers(c.ksize, read)
        for kmer in kmers:
            if kmer in self.inner_kmers: 
                self.inner_kmers[kmer] += 1

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        d = 'simulatiom' if c.simulation else 'output'
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../' + d + '/' + bed_file_name + '/' + str(c.ksize) + '/', self.previous_job_name[:-1]))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class LocalUniqueKmersGenotypingJob(map_reduce.BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = LocalUniqueKmersGenotypingJob(job_name = 'Genotyping_', previous_job_name = 'Genotyping_', batch_file_prefix = 'outer_', previous_job_batch_file_prefix = 'experimental_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        c = config.Configuration()
        likelihood = {}
        distribution = {
            '(1, 1)': (lambda x: abs(x - 2)),
            '(1, 0)': (lambda x: abs(x - 1)),
            '(0, 0)': (lambda x: abs(x - 0)),
        }
        # all kmers remaining in the track are to be considered part of the signature, no matter the number of breakpoints
        #with open(track, 'r') as track_file:
        #    break_points = json.load(track_file)
        break_points = track
        for break_point in break_points:
            if len(break_points[break_point]['right_local_unique_kmers']) == 0:
                return None
            if len(break_points[break_point]['left_local_unique_kmers']) == 0:
                return None
            likelihood[break_point] = break_points[break_point]
            likelihood[break_point]['likelihood'] = {}
            likelihood[break_point]['mean'] = {}
            inner = likelihood[break_point]['mean']['inner'] = statistics.mean(list(map(lambda kmer: break_points[break_point]['inner_kmers'][kmer], break_points[break_point]['inner_kmers'])))
            right = likelihood[break_point]['mean']['right'] = statistics.mean(list(map(lambda kmer: break_points[break_point]['right_local_unique_kmers'][kmer], break_points[break_point]['right_local_unique_kmers'])))
            left =  likelihood[break_point]['mean']['left']  = statistics.mean(list(map(lambda kmer: break_points[break_point]['left_local_unique_kmers'][kmer], break_points[break_point]['left_local_unique_kmers'])))
            #print(likelihood[break_point]['mean'])
            if inner == 0:
                likelihood[break_point]['genotype'] = '(0, 0)'
            if right == 0 or left == 0:
                likelihood[break_point]['genotype'] = '(1, 1)'
            else:
                m = statistics.mean([right, left])
                r = float(inner) / float(m)
                for zyg in distribution:
                    likelihood[break_point]['likelihood'][zyg] = distribution[zyg](r)
                choice = min(likelihood[break_point]['likelihood'].items(), key = operator.itemgetter(1))[0]
                print(choice)
                likelihood[break_point]['genotype'] = choice
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + 'genotype_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(likelihood, json_file, sort_keys = True, indent = 4)
        return path

    def reduce(self):
        c = config.Configuration()
        output = map_reduce.Job.reduce(self)
        with open(os.path.join(self.get_current_job_directory(), 'outer_merge.bed'), 'w') as bed_file:
            for track in output:
                print(output[track])
                with open(output[track], 'r') as track_file:
                    break_points = json.load(track_file)
                    chrom = track.split('_')[0]
                    begin = track.split('_')[1]
                    end = track.split('_')[2]
                    for break_point in break_points:
                        #bed_file.write('{:5}'.format(chrom) + '    ' + '{:12}'.format(begin) + '    ' + '{:12}'.format(end) + '    ' + '{:5}'.format(str(break_points[break_point]['mean']['inner'])) + '    ' +
                        #        '{:5}'.format(str(break_points[break_point]['mean']['left'])) + '    ' + '{:5}'.format(str(break_points[break_point]['mean']['right'])) + '    ' +
                        #    break_points[break_point]['genotype'] +
                        #    '\n')
                        bed_file.write(chrom + '\t' + begin + '\t' + end + '\t' + str(break_points[break_point]['mean']['inner']) + '\t' +
                            str(break_points[break_point]['mean']['left']) + '\t' + str(break_points[break_point]['mean']['right']) + '\t' +
                            break_points[break_point]['genotype'] +
                            '\n')
        return output

    def plot(self, tracks):
        left = []
        right = []
        inner = []
        for track in tracks:
            with open(tracks[track], 'r') as track_file:
                break_points = json.load(track_file)
                for break_point in break_points:
                    inner.append(break_points[break_point]['mean']['inner'])
                    right.append(len(break_points[break_point]['right_local_unique_kmers']))
                    left.append(len(break_points[break_point]['left_local_unique_kmers']))
        visualizer.histogram(inner, 'inner_kmers', self.get_current_job_directory(), x_label = 'mean of the number of reads supporting each kmer', y_label = 'number of events')
        visualizer.histogram(right, 'num_right_local_unique_kmers', self.get_current_job_directory(), x_label = 'number of unique right kmers', y_label = 'number of events')
        visualizer.histogram(left, 'num_left_local_unique_kmers', self.get_current_job_directory(), x_label = 'number of unique left kmers', y_label = 'number of events')
        #visualizer.histogram(right, 'num_right_local_unique_kmers', self.get_current_job_directory(), x_label = 'mean of the number of reads supporting each kmer', y_label = 'number of events')
        #visualizer.histogram(left, 'num_left_local_unique_kmers', self.get_current_job_directory(), x_label = 'mean of the number of reads supporting each kmer', y_label = 'number of events')
        self.plot_accuracy_over_event_length(tracks)


    def plot_accuracy_over_event_length(self, tracks):
        bins = {}
        L = []
        for track in tracks:
            print(track)
            tokens = track.split('_')
            start = int(tokens[1])
            end = int(tokens[2])
            l = end - start
            L.append(l)
            with open(tracks[track], 'r') as track_file:
                break_points = json.load(track_file)
                for break_point in break_points:
                    b = 100 * int(l / 100)
                    if not b in bins:
                        bins[b] = []
                    if break_points[break_point]['genotype'] == '(0, 0)':
                        bins[b].append(0)
                    else:
                        bins[b].append(1)
        bins = sorted(bins.items(), key = operator.itemgetter(0))
        visualizer.bar(x = list(map(lambda b: b[0], bins)), ys = [list(map(lambda b: len(b[1]), bins)), list(map(lambda b: sum(b[1]), bins))], path = self.get_current_job_directory(), name = 'accuracy_per_event_length', x_label = 'event lenght rounded to the nearest hundred', y_label = 'genotypinh accuracy')
        visualizer.histogram(x = L, name = 'event_length_distribution', path = self.get_current_job_directory(), x_label = 'event length', y_label = 'number of events')


# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenotypingAnalysisJob(map_reduce.BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = GenotypingAnalysisJob(job_name = 'Genotyping_', previous_job_name = 'Genotyping_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = 1

    def load_inputs(self):
        self.true_positive = {}
        self.true_negative = {}
        self.false_positive = {}
        self.false_negative = {}
        classes = ['true_positive', 'true_negative', 'false_positive', 'false_negative']
        for c in classes:
            setattr(self, c, {})
            with open(os.path.join(self.get_current_job_directory(), c + '.bed'), 'r') as bed_file:
                for track in bed.parse_bed_file(bed_file):
                    track_name = track[0] + '_' + track[1] + '_' + track[2]
                    getattr(self, c)[track_name] = track
        for c in classes:
            self.plot_kmer_count_variation(c, getattr(self, c))
        exit()

    def plot_kmer_count_variation(self, batch_name, batch):
        variance = {'inner': [], 'novel': []}
        for track_name in batch:
            path = os.path.join(self.get_current_job_directory(), 'genotype_' + track_name + '.json')
            if not os.path.exists(path):
                print(red('not found:', path))
                continue
            with open(path, 'r') as json_file:
                break_points = json.load(json_file)
                for break_point in break_points:
                    inner_kmers = list(map(lambda x: x[1], break_points[break_point]['inner_kmers'].items()))
                    novel_kmers = list(map(lambda x: x[1], break_points[break_point]['novel_kmers'].items()))
                    if len(inner_kmers):
                        variance['inner'].append(statistics.variance(inner_kmers))
                    if len(novel_kmers):
                        variance['novel'].append(statistics.variance(novel_kmers))
        v = statistics.variance(variance['inner'])
        m = statistics.mean(variance['inner'])
        data = [graph_objs.Histogram(x = variance['inner'], xbins = dict(start = 0, size = 5, end = 1000), text = 'Mean:' + str(m) + ' STD: ' + str(math.sqrt(v)))]
        plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), batch_name + '_inner.html'), auto_open = False)
        v = statistics.variance(variance['novel'])
        m = statistics.mean(variance['novel'])
        data = [graph_objs.Histogram(x = variance['novel'], xbins = dict(start = 5, size = 5, end = 1000), text = 'Mean:' + str(m) + ' STD: ' + str(math.sqrt(v)))]
        plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), batch_name + '_novel.html'), auto_open = False)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class DepthOfCoverageEstimationJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = DepthOfCoverageEstimationJob(job_name = 'DepthOfCoverageEstimationJob_', previous_job_name = "", **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare(self):
        self.kmers = {}

    def load_inputs(self):
        c = config.Configuration()
        self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        # check if we already have the kmers
        if os.path.isfile(os.path.join(self.get_current_job_directory(), 'kmers.json')):
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'r') as kmers_file:
                print('Exon kmers already extracted, reloading from cache')
                self.kmers = json.load(kmers_file)
                # setting num_threads to 0 will bypass all execution and jump to reduce
                self.num_threads = 0
                return
        # load exonic regions
        tracks = bed.read_tracks(c.exons)
        for index in range(0, self.num_threads):
            self.batch[index] = {}
        # split exonic regions into several batches
        index = 0
        for track in tracks:
            print(track)
            self.batch[index][track] = tracks[track]
            index += 1
            if index == self.num_threads:
                index = 0

    def transform(self, track, track_name):
        c = config.Configuration()
        self.kmers[track_name] = {}
        seq = bed.extract_sequence(track)
        for kmer in extract_kmers(c.ksize, seq):
            if not kmer in self.kmers[track_name]:
                count = self.counts_provider.get_kmer_count(kmer)
                self.kmers[track_name][kmer] = count
        return True

    def output_batch(self, batch):
        path = os.path.join(self.get_current_job_directory(), 'kmers_' + str(self.index) + '.json')
        with open(path, 'w') as json_file:
            json.dump(self.kmers, json_file, sort_keys = True, indent = 4)

    def reduce(self):
        if self.num_threads != 0:
            # merge all the kmer counts from previous steps
            self.kmers = {}
            for i in range(0, self.num_threads):
                print('adding batch', i)
                path = os.path.join(self.get_current_job_directory(), 'kmers_' + str(i) + '.json')
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    for track in batch:
                        for kmer in batch[track]:
                            self.kmers[kmer] = batch[track][kmer]
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
        else:
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
                self.kmers = json.load(json_file)
        #
        #print('got', len(self.kmers), 'kmers, sampling 100000 randomly')
        #sample = random.sample(self.kmers.items(), min(len(self.kmers), 100000))
        #for kmer in self.kmers:
        #    self.kmers[kmer] = self.get_kmer_count(kmer, self.index, False)
        # calculate mean and std
        self.counts = list(map(lambda x: self.kmers[x], list(self.kmers.keys())))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print('mean:', self.mean)
        print('std:', self.std)
        # filter outliers
        self.counts = list(filter(lambda x: x < 3 * self.mean, self.counts))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print('mean:', self.mean)
        print('std:', self.std)
        # filter anything appearing more than twice the medium, 4x coverage or more, repeatimg kmer
        self.counts = list(filter(lambda x: x < 2 * self.mean, self.counts))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print('mean:', self.mean)
        print('std:', self.std)
        #
        with open(os.path.join(self.get_current_job_directory(), 'stats.json'), 'w') as json_file:
            json.dump({ 'mean': self.mean, 'std': self.std }, json_file, sort_keys = True, indent = 4)

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_current_job_directory(self):
        c = config.Configuration()
        if c.simulation:
            return map_reduce.Job.get_current_job_directory(self)
        else:
            return self.get_output_directory()

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
