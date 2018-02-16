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

from kmer.sv import StructuralVariation, Inversion, Deletion
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
# MapReduce job for finding StructuralVariation breakpoints
# Algorithm: starts with a set of structural variation events and their approximate breakpoints and tries to refine them
# considers a radius of [-50, 50] around each end of the breakpoint, and for each pair of endpoints within that radius considers
# the area as the structural variations and applies it to the reference genome to generate a set of kmers. Discards those endpoints
# whose kmers do not all appear in the base genome the event was detected in. 
# Output: Reports the remaining boundary candidates with their list of associated kmers and the count of those kmers in the
# base genome.
# ============================================================================================================================ #
# ============================================================================================================================ #

class BreakPointJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = BreakPointJob(job_name = 'break_point_', previous_job_name = '', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed the BED file from which to read the tracks (default: CHM1.Lumpy)
        # --threads the number of threads to use (default 48)
        # --reference the reference assemby on which the coordinates in the BED are based (default: hg19)
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.bedtools = pybedtools.BedTool(c.bed_file)
        self.num_threads = max(c.max_threads, len(self.bedtools))

    def load_inputs(self):
        c = config.Configuration()
        self.radius = 50
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

    def run_batch(self, batch):
        c = config.Configuration()
        sv_type = self.get_sv_type()
        output = {}
        n = 0
        start = time.time()
        for track in batch:
            name = re.sub(r'\s+', '_', str(track).strip()).strip()
            # track coordinates might exceed boundaries of the chromosome
            try:
                sv = sv_type(track = track, radius = self.radius)
                output[name] = self.transform(sv)
            except pybedtools.helpers.BEDToolsError as e:
                print(e)
            t = time.time()
            c = float(n) / len(batch)
            print('index:', self.index, 'completion:', c, 'ETA:', ((1.0 - c) * (t - start) / c) / 3600, 'hours')
        self.output_batch(output)
        print(colorama.Fore.GREEN + 'process ', self.index, ' done')

    def transform(self, sv):
        c = config.Configuration()
        break_points = self.extract_boundary_kmers(sv)
        break_points = self.calc_break_point_scores(break_points, sv)
        for break_point in break_points:
            # counts for reference not available at this time, need a differnt count_sever instance, need to calc separately
            #for kmer in break_points[break_point]['reference_kmers']:
            #    break_points[break_point]['reference_kmers'][kmer] = -1
            for kmer in break_points[break_point]['kmers']:
                break_points[break_point]['kmers'][kmer] = count_server.get_kmer_count(kmer, self.index, False)
            # save the number of boundary candidates
        return break_points

    def reduce(self):
        c = config.Configuration()
        output = {}
        path = os.path.join(self.get_current_job_directory(), 'break_points')
        if not os.path.exists(path):
            os.makedirs(path)
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    for track in batch:
                        with open(os.path.join(self.get_current_job_directory(), 'break_points', track + '.json'), 'w') as track_file:
                            json.dump({track: batch[track]}, track_file, sort_keys = True, indent = 4)
                    output.update(batch)
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)

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
        for begin in range(-self.radius, self.radius + 1) :
            for end in range(-self.radius, self.radius + 1) :
                kmers, boundary = sv.get_signature_kmers(begin, end)
                if not kmers:
                    # skip this candidate, has overlapping ends
                    continue
                #reference_kmers = sv.get_reference_signature_kmers(begin, end)
                break_points['(' + str(begin) + ',' + str(end) + ')'] = {
                    'boundary': boundary,
                    'kmers': kmers,
                    #'reference_kmers': reference_kmers
                }
        return break_points

    # prunes a break points if not all its kmers appear in the counttable
    def calc_break_point_scores(self, break_points, sv):
        c = config.Configuration()
        for break_point in break_points:
            n = 0
            for kmer in break_points[break_point]['kmers']:
                count = count_server.get_kmer_count(kmer, self.index, False)
                if count != 0:
                    n = n + 1
            break_points[break_point]['score'] = float(n) / float(len(break_points[break_point]['kmers']))
        return break_points

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

class AllBreakPointsJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = AllBreakPointsJob(job_name = 'AllBreakPoints_', previous_job_name = '')
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self):
        # --bed the BED file from which to read the tracks (sensible default provided)
        # --threads the number of threads to use (sensible default provided)
        # --reference the reference assemby on which the coordinates in the BED are based
        pass

    def find_thread_count(self):
        pass

    def load_inputs(self):
        c = config.Configuration()
        bedtools = pybedtools.BedTool(c.bed_file)
        self.radius = 50
        # split variations into batches
        n = 0
        for track in bedtools:
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
        self.num_threads = len(self.batch)
        print('running on ', self.num_threads, ' threads')

    def run_batch(self, batch):
        c = config.Configuration()
        sv_type = self.get_sv_type()
        output = {}
        n = 0
        start = time.time()
        for track in batch:
            name = re.sub(r'\s+', '_', str(track).strip()).strip()
            # track coordinates might exceed boundaries of the chromosome
            try:
                sv = sv_type(track = track, radius = self.radius)
                output[name] = self.transform(sv)
            except pybedtools.helpers.BEDToolsError as e:
                print(e)
            n = n + 1
            t = time.time()
            p = float(n) / len(batch)
            eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
            print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
        self.output_batch(output)
        print(colorama.Fore.GREEN + 'process ', self.index, ' done')

    def transform(self, sv):
        c = config.Configuration()
        break_points = self.generate_break_points(sv)
        return break_points

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

    def generate_break_points(self, sv):
        c = config.Configuration()
        break_points = {}
        for begin in range(-self.radius, self.radius + 1) :
            for end in range(-self.radius, self.radius + 1) :
                kmers, boundary = sv.get_signature_kmers(begin, end)
                if not kmers:
                    # skip this candidate, has overlapping ends
                    continue
                #reference_kmers = sv.get_reference_signature_kmers(begin, end)
                break_points['(' + str(begin) + ',' + str(end) + ')'] = {
                    'boundary': boundary,
                    'kmers': kmers,
                    #'reference_kmers': reference_kmers
                }
        return break_points

# ============================================================================================================================ #
# ============================================================================================================================ #
# Counts the kmers within the boundaries of each structural variation exactly. 
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExactCountBreakPointsJob(map_reduce.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        for i in range(0, 24):
            job = ExactCountBreakPointsJob(job_name = 'ExactCountBreakPoints_', previous_job_name = 'AllBreakPoints_', input_chunk = i, input_chunk_size = 2)
            job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        tracks = {}
        for batch in range(self.input_chunk * self.input_chunk_size, (self.input_chunk + 1) * self.input_chunk_size):
            print('loading batch', batch)
            path = os.path.join(self.get_previous_job_directory(), 'batch_' + str(batch) + '.json')
            with open(path, 'r') as json_file:
                tracks.update(json.load(json_file))
        # accumulate all kmers from all tracks into a large dict to improve performance
        # if we keep each tracks's kmers in a seperate dict, runtime will grow linearly with number of tracks
        self.kmers = {}
        for track in tracks:
            # remember these kmers already come in the canonical representation
            for break_point in tracks[track]:
                for kmer in tracks[track][break_point]['kmers']:
                    self.kmers[kmer] = 0
        print('done loading kmers,', len(self.kmers), 'total')
        # these are just dummies
        for index in range(0, self.num_threads):
            self.batch[index] = {} # avoid overrding extra methods from MapReduce

    def reduce(self):
        kmers = self.merge_counts()
        # reload tracks
        path = os.path.join(self.get_current_job_directory(), 'merge_' + str(self.input_chunk) + '.json')
        with open(path, 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)
        """
        racks = json.load(json_file)
        # update exact counts for break points
        for track in tracks:
            for break_point in tracks[track]:
                for kmer in tracks[track][break_point]['kmers']:
                    tracks[track][break_point]['kmers'] = kmers[kmer]
        # prune low quality break points
        for track in tracks:
            for break_point in tracks[track]:
                self.calc_break_point_score(break_point)
        """

    def calc_break_point_score(self, break_point):
        c = config.Configuration()
        n = 0
        for kmer in break_point['kmers']:
            if break_point['kmers'][kmer] != 0:
                n = n + 1
        break_point['score'] = float(n) / float(len(break_point['kmers']))

# ============================================================================================================================ #
# ============================================================================================================================ #
# Utitlizes the likelihood model to find the most likely breakpoint for each structural variation
# ============================================================================================================================ #
# ============================================================================================================================ #

class MostLikelyBreakPointsJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = MostLikelyBreakPointsJob(job_name = 'MostLikelyBreakPoints_', previous_job_name = 'CountKmersExactJob_', **kwargs)
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
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            self.tracks = json.load(json_file)
        for index in range(0, self.num_threads):
            self.batch[index] = {}
        index = 0
        for track in self.tracks:
            self.batch[index][track] = self.tracks[track]
            index += 1
            if index == self.num_threads:
                index = 0

    def transform(self, track, track_name):
        c = config.Configuration()
        self.coverage = c.coverage
        self.std = c.std
        distribution = {
            '(1, 1)': statistics.NormalDistribution(mean = self.coverage, std = self.std),
            '(1, 0)': statistics.NormalDistribution(mean = self.coverage / 2, std = self.std),
            '(0, 0)': statistics.ErrorDistribution(1.0 / 1000)
        }
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
        # calculate likelihoods
        novel_kmers = track['novel_kmers']
        for kmer in novel_kmers:
            for break_point in likelihood['break_points']:
                r = 0
                if break_point in novel_kmers[kmer]['break_points']:
                    r = distribution['(1, 1)'].log_pmf(novel_kmers[kmer]['actual_count'])
                    likelihood['break_points'][break_point]['novel_kmers'][kmer] = novel_kmers[kmer]#{'count': novel_kmers[kmer]['actual_count'], 'likelihood': r}
                else:
                    r = distribution['(0, 0)'].log_pmf(novel_kmers[kmer]['actual_count'])
                likelihood['break_points'][break_point]['likelihood'] += r
        # find the maximum one and keep it easily accessible in the output
        # TODO sort and pic 5 best
        l = list(map(lambda x: (x, likelihood['break_points'][x]['likelihood']), likelihood['break_points']))
        m = max(l, key = operator.itemgetter(1))[0]
        likelihood['most_likely'] = {m: likelihood['break_points'][m]}
        return likelihood

    def reduce(self):
        c = config.Configuration()
        output = {}
        path = (os.path.join(self.get_current_job_directory(), 'most_likely'))
        if not os.path.exists(path):
            os.makedirs(path)
        distribution = []
        # collect all most likely breakpoints in a single for reference
        most_likely = {}
        # dump each event's most likely break point to a separate file for easier access
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    for track in batch:
                        most_likely[track] = {}
                        most_likely[track]['break_points'] = {}
                        path = os.path.join(self.get_current_job_directory(), 'most_likely', track + '.most_likely.json')
                        with open(path, 'w') as most_likely_file:
                            json.dump(batch[track]['most_likely'], most_likely_file, sort_keys = True, indent = 4)
                            # this prints extra stuff for debugging purposes
                            for break_point in batch[track]['most_likely']:
                                most_likely[track]['break_points'][break_point] = batch[track]['most_likely'][break_point]
                                break
                                novel_kmers = batch[track]['most_likely'][break_point]['novel_kmers']
                                n = 0
                                m = 0
                                k = 0
                                for kmer in novel_kmers:
                                    if abs(novel_kmers[kmer]['count'] - c.coverage) <= c.std:
                                        n += 1
                                    if abs(novel_kmers[kmer]['count'] - c.coverage) <= 2 * c.std:
                                        m += 1
                                    if abs(novel_kmers[kmer]['count'] - c.coverage) <= 3 * c.std:
                                        k += 1
                                distribution.append(float(n) / len(novel_kmers))
                                #print('{:30}'.format(track), '{:12}'.format(break_point), '{:4d}'.format(len(novel_kmers)), '{:4d}'.format(n), '{:4d}'.format(m), '{:4d}'.format(k))
                    output.update(batch)
            print('merging,', i, ' out of', self.num_threads, ' done')
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'most_likely.json'), 'w') as json_file:
            json.dump(most_likely, json_file, sort_keys = True, indent = 4)
        self.plot(distribution)

    def plot(self, x):
        path = os.path.join(self.get_current_job_directory(), 'plot_most_likely_break_point_kmer_counts.html')
        trace = graph_objs.Histogram(
            x = x,
            histnorm = 'count',
            xbins = dict(
                start = 0.0,
                end = 1.0,
                size = 0.05
            )
        )
        layout = graph_objs.Layout(title = 'Deviation from Mean')
        fig = graph_objs.Figure(data = [trace], layout = layout)
        plotly.plot(fig, filename = path)

# ============================================================================================================================ #
# Helper job for rapidly generating likelihood plots
# ============================================================================================================================ #

class MostLikelyBreakPointsPlottingJob(map_reduce.Job):

    @staticmethod
    def launch():
        job = MostLikelyBreakPointsPlottingJob(job_name = 'MostLikelyBreakPoints_', previous_job_name = 'MostLikelyBreakPoints_')
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
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    #
    #ExactCountBreakPointsJob.launch()
    #AllBreakPointsJob.launch()
    #MostLikelyBreakPointsJob.launch(resume_from_reduce = True)
    #MostLikelyBreakPointsPlottingJob.launch()
    BreakPointJob.launch(resume_from_reduce = True)
