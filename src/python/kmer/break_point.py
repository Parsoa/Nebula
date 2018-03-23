import io
import os
import re
import pwd
import sys
import copy
#import json
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
    count_server,
)

from kmer.sv import StructuralVariation, Inversion, Deletion, SNP
from kmer.kmers import *
from kmer.commons import *
from kmer.commons import pretty_print as print

import khmer
import colorama
import pybedtools

import rapidjson as json
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
        job = BreakPointJob(job_name = 'ExtractBreakPointsJob_', previous_job_name = '', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed the BED file from which to read the tracks (default: CHM1.Lumpy)
        # --snp the BED file containing the coordinates of commons SNPs (optional)
        # --threads the number of threads to use (default 48)
        # --reference the reference assemby on which the coordinates in the BED are based (default: hg19)
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        c = config.Configuration()
        self.bedtools = pybedtools.BedTool(c.bed_file)
        self.radius = 50
        self.snps = {}
        # SNPs are indexed by their position
        for track in bed.parse_bed_file(open(c.snp, 'r')):
            self.snps[track[2]] = SNP(chrom = track[1], begin = track[2], end = track[3], variants = track[9])
        # split events into batches
        n = 0
        for track in self.bedtools:
            name = re.sub(r'\s+', '_', str(track).strip()).strip()
            # too large, skip
            if track.end - track.start > 1000000:
                print(colorama.Fore.RED + 'skipping ', name, ', too large')
                continue
            index = n % c.max_threads 
            if not index in self.batch:
                self.batch[index] = []
            self.batch[index].append(track)
            print(colorama.Fore.BLUE + 'assigned ', name, ' to ', index)
            n = n + 1
            self.num_threads = max(self.num_threads, index + 1)

    def run_batch(self, batch):
        c = config.Configuration()
        sv_type = self.get_sv_type()
        self.tracks = {}
        n = 0
        start = time.time()
        for track in batch:
            name = re.sub(r'\s+', '_', str(track).strip()).strip()
            # track coordinates might exceed boundaries of the chromosome
            try:
                sv = sv_type(track = track, radius = self.radius)
                break_points = self.transform(sv)
                self.output_break_points(name, break_points)
            except pybedtools.helpers.BEDToolsError as e:
                print(e)
            n = n + 1
            t = time.time()
            c = float(n) / len(batch)
            print('index:', self.index, 'completion:', c, 'ETA:', ((1.0 - c) * (t - start) / c) / 3600, 'hours')
        self.output_batch(self.tracks)

    def transform(self, sv):
        c = config.Configuration()
        break_points = self.extract_boundary_kmers(sv)
        if not break_points:
            return None
        for break_point in break_points:
            for kmer in break_points[break_point]['kmers']:
                break_points[break_point]['kmers'][kmer] = count_server.get_kmer_count(kmer, self.index, False)
            self.calc_break_point_score(break_points[break_point])
        return break_points

    def reduce(self):
        c = config.Configuration()
        output = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    output.update(batch)
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)
        self.sort(output)
        self.plot(output)

    def output_break_points(self, track_name, break_points):
        path = os.path.join(self.get_current_job_directory(), 'break_points_' + track_name  + '.json') 
        self.tracks[track_name] = path
        json_file = open(path, 'w')
        json.dump({'break_points': break_points}, json_file, sort_keys = True, indent = 4)
        json_file.close()

    # ============================================================================================================================ #
    # job-specific helpers
    # ============================================================================================================================ #

    def get_sv_type(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        sv_type = bed_file_name.split('.')[-2]
        if sv_type == 'DEL':
            return Deletion
        if sv_type == 'INV':
            return Inversion
        return StructuralVariation

    def extract_boundary_kmers(self, sv):
        c = config.Configuration()
        break_points = {}
        if c.snp:
            events = sv.find_snps_within_boundaries(self.snps)
            print(len(events), 'events found')
        for begin in range(-self.radius, self.radius + 1) :
            for end in range(-self.radius, self.radius + 1) :
                kmers, boundary = sv.get_signature_kmers(begin, end)
                if not kmers:
                    continue
                break_points['(' + str(begin) + ',' + str(end) + ')'] = {
                    'boundary': boundary,
                    'kmers': kmers,
                }
                if c.snp:
                    for event in events:
                        for kmers, boundary, variant in sv.get_signature_kmers_with_variation(events[event], begin, end):
                            break_points['(' + str(begin) + ',' + str(end) + ')_' + events[event].begin + '_' + variant] = {
                                'boundary': boundary,
                                'kmers': kmers,
                            }
        print(len(break_points), 'break points found')
        return break_points if break_points else None

    # prunes a break points if not all its kmers appear in the counttable
    def calc_break_point_score(self, break_point):
        c = config.Configuration()
        n = len(list(filter(lambda kmer: break_point['kmers'][kmer] != 0, break_point['kmers'])))
        break_point['score'] = float(n) / float(len(break_point['kmers']))

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
        job = MostLikelyBreakPointsJob(job_name = 'MostLikelyBreakPointsJob_', previous_job_name = 'CountKmersExactJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed the BED file including te current set of structural variations
        # --std the standard deviation of the genome these structural variations were extracted from
        # --threads the number of processes to spawn
        # --coverage the average depth of coevrage the genome these structural variations were extracted from
        pass

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
        self.coverage = c.coverage
        self.std = c.std
        self.distribution = {
            '(1, 1)': statistics.NormalDistribution(mean = self.coverage, std = self.std),
            '(1, 0)': statistics.NormalDistribution(mean = self.coverage / 2, std = self.std),
            '(0, 0)': statistics.ErrorDistribution(1.0 / 1000)
        }

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(track, 'r') as track_file:
            track = json.load(track_file)
        # find all the break points
        likelihood = {
            'break_points': {}
        }
        for kmer in track['novel_kmers']:
            for break_point in track['novel_kmers'][kmer]['break_points']:
                if not break_point in likelihood['break_points']:
                    likelihood['break_points'][break_point] = {
                        'likelihood': 0,
                        'novel_kmers': {} # we are also interested in knowing the kmer with for each one
                    }
        # in case no breakpoints with adequate score existed for this track
        if len(likelihood['break_points']) == 0:
            return None
        print(track_name, ': ', len(likelihood['break_points']), ' break points')
        # calculate likelihoods
        m = None
        novel_kmers = track['novel_kmers']
        for kmer in novel_kmers:
            r_1_1 = self.distribution['(1, 1)'].log_pmf(novel_kmers[kmer]['actual_count'])
            r_0_0 = self.distribution['(0, 0)'].log_pmf(novel_kmers[kmer]['actual_count'])
            for break_point in likelihood['break_points']:
                if break_point in novel_kmers[kmer]['break_points']:
                    likelihood['break_points'][break_point]['novel_kmers'][kmer] = {
                        'count': novel_kmers[kmer]['count'],
                        'likelihood': r_1_1,
                        'actual_count': novel_kmers[kmer]['actual_count']
                    }
                    likelihood['break_points'][break_point]['likelihood'] += r_1_1
                else:
                    likelihood['break_points'][break_point]['likelihood'] += r_0_0
                m = break_point if not m or likelihood['break_points'][break_point]['likelihood'] > likelihood['break_points'][m]['likelihood'] else m
        # find the maximum one and keep it easily accessible in the output
        # l = list(map(lambda x: (x, likelihood['break_points'][x]['likelihood']), likelihood['break_points']))
        # m = max(l, key = operator.itemgetter(1))[0]
        path = os.path.join(self.get_current_job_directory(), 'likelihood_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(likelihood, json_file, sort_keys = True, indent = 4)
        path = os.path.join(self.get_current_job_directory(), 'most_likely_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({m: likelihood['break_points'][m]}, json_file, sort_keys = True, indent = 4)
        return path

    def plot(self, x):
        return
        path = os.path.join(self.get_current_job_directory(), 'plot_most_likely_break_point_kmer_counts.html')
        trace = graph_objs.Histogram(
            x = x,
            histnorm = 'count',
            xbins = dict(
                start = 0.0,
                end = 1.0,
                sizs = 0.05
            )
        )
        layout = graph_objs.Layout(title = 'Deviation from Mean')
        fig = graph_objs.Figure(data = [trace], layout = layout)
        plotly.plot(fig, filename = path)

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
        print('threads:', colorama.Fore.GREEN, self.num_threads)

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
        print(colorama.Fore.CYAN + 'reduce')
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
        BreakPointJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'MostLikelyBreakPointsJob':
        MostLikelyBreakPointsJob.launch(resume_from_reduce = c.resume_from_reduce)
