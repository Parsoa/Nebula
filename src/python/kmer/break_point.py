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

from kmer import (
    bed,
    sets,
    config,
    counttable,
    map_reduce,
    statistics,
)

from kmer.sv import StructuralVariation, Inversion, Deletion, SNP
from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import pybedtools

import plotly.offline as plotly
import plotly.graph_objs as graph_objs

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

class ExtractBreakPointsJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExtractBreakPointsJob(job_name = 'ExtractBreakPointsJob_', previous_job_name = '', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        c = config.Configuration()
        self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.bedtools = pybedtools.BedTool(c.bed_file)
        round_robin(self.bedtools, lambda track: re.sub(r'\s+', '_', str(track).strip()).strip(), lambda track: track.end - track.start > 1000000) 

    def transform(self, track, track_name):
        sv = sv_type(track = track)
        c = config.Configuration()
        break_points = self.extract_break_points(sv, track_name)
        if not break_points:
            return None
        path = os.path.join(self.get_current_job_directory(), 'break_points_' + track_name  + '.json') 
        json_file = open(path, 'w')
        json.dump({'break_points': break_points}, json_file, sort_keys = True, indent = 4)
        json_file.close()
        return path

    # ============================================================================================================================ #
    # job-specific helpers
    # ============================================================================================================================ #

    def get_sv_type(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        if bed_file_name.find('DEL') != -1:
            return Deletion
        if bed_file_name.find('INV') != -1:
            return Inversion
        return Deletion

    def extract_break_points(self, sv):
        c = config.Configuration()
        break_points = {}
        inner_kmers = {}
        local_unique_kmers = {}
        for kmer in sv.get_inner_kmers(self.reference_counts_provider.get_kmer_count):
            inner_kmers[kmer] = self.counts_provider.get_kmer_count(kmer)
        if len(inner_kmers) == 0:
            local_unique_kmers = sv.get_local_unique_kmers(self.reference_counts_provider.get_kmer_count, c = 5, n = 1000)
            if len(local_unique_kmers) == 0:
                print(red(tr
            for kmer in sv.get_near_boundary_inner_kmers():
                inner_kmers[kmer] = self.counts_provider.get_kmer_count(kmer)
        for begin in range(-c.radius, c.radius + 1):
            for end in range(-c.radius, c.radius + 1):
                boundary_kmers, boundary = sv.get_signature_kmers(begin, end)
                if len(boundary_kmers) == 0:
                    continue
                name = '(' + str(begin) + ',' + str(end) + ')'
                break_points[name] = {
                    'boundary': boundary,
                    'inner_kmers': inner_kmers,
                    'local_unique_kmers': local_unique_kmers
                }
                novel_kmers = list(filter(lambda kmer: self.reference_counts_provider.get_kmer_count(kmer) == 0, boundary_kmers)) 
                break_points[name]['novel_kmers'] = {kmer: self.counts_provider.get_kmer_count(kmer) for kmer in novel_kmers}
                score = float(len(break_points[name]['novel_kmers'])) / len(boundary_kmers)
                break_points[name]['score'] = score
        return break_points if break_points else None

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Utitlizes the likelihood model to find the most likely breakpoint for each structural variation
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class MostLikelyBreakPointsJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = MostLikelyBreakPointsJob(job_name = 'MostLikelyBreakPointsJob_', previous_job_name = 'ExtractBreakPointsJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        c = config.Configuration()
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            paths = json.load(json_file)
        for index in range(0, self.num_threads):
            self.batch[index] = {}
        index = 0
        for track in paths:
            self.batch[index][track] = paths[track]
            index += 1
            if index == self.num_threads:
                index = 0
        self.distribution = {
            '(1, 1)': statistics.NormalDistribution(mean = c.coverage, std = c.std),
            '(1, 0)': statistics.NormalDistribution(mean = c.coverage / 2, std = c.std),
            '(0, 0)': statistics.ErrorDistribution(1.0 / 1000),
        }

    def transform(self, track, track_name):
        c = config.Configuration()
        print(track_name)
        with open(track, 'r') as track_file:
            break_points = json.load(track_file)['break_points']
        # find all the break points
        likelihood = {
            'break_points': {}
        }
        # accumulate all kmers in one dictionary
        kmers = {}
        for break_point in break_points:
            for kmer in break_points[break_point]['novel_kmers']:
                count = break_points[break_point]['novel_kmers'][kmer]
                if count:
                    kmers[kmer] = count
        # calculate likelihoods
        m = None
        for break_point in break_points:
            likelihood['break_points'][break_point] = {}
            likelihood['break_points'][break_point]['likelihood'] = 0
            likelihood['break_points'][break_point]['novel_kmers'] = break_points[break_point]['novel_kmers']
            likelihood['break_points'][break_point]['inner_kmers'] = break_points[break_point]['inner_kmers']
            for kmer in kmers:
                count = kmers[kmer]
                r_1_1 = self.distribution['(1, 1)'].log_pmf(count)
                r_1_0 = self.distribution['(1, 0)'].log_pmf(count)
                r_0_0 = self.distribution['(0, 0)'].log_pmf(count)
                if kmer in break_points[break_point]['novel_kmers']:
                    likelihood['break_points'][break_point]['likelihood'] += r_1_1
                else:
                    likelihood['break_points'][break_point]['likelihood'] += r_0_0
            m = break_point if not m or likelihood['break_points'][break_point]['likelihood'] > likelihood['break_points'][m]['likelihood'] else m
        if not m:
            print(red(track_name), 'no novel kmers found for track', red(track_name))
            return None
        path = os.path.join(self.get_current_job_directory(), 'likelihood_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(likelihood, json_file, sort_keys = True, indent = 4)
        path = os.path.join(self.get_current_job_directory(), 'most_likely_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({m: likelihood['break_points'][m]}, json_file, sort_keys = True, indent = 4)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# Helper job for rapidly generating likelihood plots
# ============================================================================================================================ #
# ============================================================================================================================ #

class MostLikelyBreakPointsPlottingJob(map_reduce.Job):

    @staticmethod
    def launch():
        job = MostLikelyBreakPointsPlottingJob(job_name = 'MostLikelyBreakPointsJob_', previous_job_name = 'MostLikelyBreakPointsJob_')
        job.execute()

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads
        print('threads:', green(self.num_threads))

    def load_inputs(self):
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            self.tracks = json.load(json_file)
        for index in range(0, self.num_threads):
            self.batch[index] = {}
        index = 0
        print(len(self.tracks))
        for track in self.tracks:
            self.batch[index][track] = self.tracks[track]
            index += 1
            if index == self.num_threads:
                index = 0

    def transform(self, track, track_name):
        self.plot_sorted_break_point_likelihoods(track, track_name)
        self.plot_likelihood_heatmap(track, track_name)
        self.plot_kmer_count_heatmap(track, track_name)
        return None

    def output_batch(self, batch):
        # no need to output anything
        pass

    def plot_kmer_count_heatmap(self, track, track_name):
        self.radius = 50
        x = []
        for begin in range(-self.radius, self.radius + 1) :
            x.append([])
            for end in range(-self.radius, self.radius + 1) :
                break_point = '(' + str(begin) + ',' + str(end) + ')'
                if break_point in track['break_points']:
                    x[begin + self.radius].append(len(track['break_points'][break_point]['novel_kmers']))
                else:
                    x[begin + self.radius].append(0)
        path = os.path.join(self.get_current_job_directory(), track_name + '_novel_kmer_count.html')
        trace = graph_objs.Heatmap(z = x)
        data = [trace]
        plotly.plot(data, filename = path, auto_open = False)

    def plot_likelihood_heatmap(self, track, track_name):
        self.radius = 50
        z = []
        x = []
        y = []
        m = max(list(map(lambda x: track['break_points'][x]['likelihood'], track['break_points'])))
        for begin in range(-self.radius, self.radius + 1) :
            z.append([])
            y.append(begin)
            x = []
            for end in range(-self.radius, self.radius + 1) :
                x.append(end)
                break_point = '(' + str(begin) + ',' + str(end) + ')'
                if break_point in track['break_points']:
                    z[begin + self.radius].append(track['break_points'][break_point]['likelihood'])
                else:
                    z[begin + self.radius].append(m + 1000)
        path = os.path.join(self.get_current_job_directory(), track_name + '_likelihood_heatmap.html')
        trace = graph_objs.Heatmap(z = z, x = x, y = y)
        data = [trace]
        plotly.plot(data, filename = path, auto_open = False)

    def plot_sorted_break_point_likelihoods(self, track, track_name):
        likelihoods = sorted(list(map(lambda x: track['break_points'][x]['likelihood'], track['break_points'])))
        x = list(range(0, len(likelihoods)))
        path = os.path.join(self.get_current_job_directory(), track_name + '_likelihood_sorted.html')
        trace = graph_objs.Scatter(x = x, y = likelihoods, mode = 'lines')
        data = [trace]
        plotly.plot(data, filename = path, auto_open = False)

    def reduce(self):
        exit()
        # no need to merge anything

    def get_current_job_directory(self):
        return os.path.abspath(os.path.join(self.get_output_directory(), self.job_name[:-1], 'plots'))

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
    if c.job == 'ExtractBreakPointsJob':
        ExtractBreakPointsJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'MostLikelyBreakPointsJob':
        MostLikelyBreakPointsJob.launch(resume_from_reduce = c.resume_from_reduce)
