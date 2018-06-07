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
import random
import argparse
import traceback

from kmer import (
    bed,
    sets,
    config,
    counttable,
    map_reduce,
    statistics,
)
from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import numpy
#import rapidjson as json
#import plotly.offline as plotly
#import plotly.graph_objs as graph_objs

# ============================================================================================================================ #
# ============================================================================================================================ #
# Once we have extracted the kmers for every breakpoint of a structural variation, we will be interested in knowing
# which of thses kmers are actually novel ones that don't normally occur in the genome.
# REQUIRES a count_server for the reference genome
# ============================================================================================================================ #
# ============================================================================================================================ #

class NovelKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = NovelKmersJob(job_name = 'NovelKmersJob_', previous_job_name = 'ExtractBreakPointsJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed to indicate the set of structural variations we are running this for (default value supplied)
        # --threads to indicate number of threads (default value provided)
        pass

    # the output from previous stage is too large, we need to load it sequentially to avoid memory allocation errors
    # json libraries don't seem to support reading sequenctially from a file so we will read each batch of the previous job
    # and export each of its tracks to a separate file in parallel.
    def load_inputs(self):
        self.batch = {}
        self.counts_provider = counttable.JellyfishCountsProvider()
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            tracks = json.load(json_file)
            for index in range(0, self.num_threads):
                self.batch[index] = {}
            index = 0
            for track in tracks:
                self.batch[index][track] = tracks[track]
                index = index + 1
                if index == self.num_threads:
                    index = 0

    def transform(self, path, track_name):
        # load each job's input separately to reduce memory footprint
        with open(path, 'r') as json_file:
            track = json.load(json_file)['break_points']
        break_points = {}
        n = 0
        for break_point in track:
            # ignore weak break_points
            if float(track[break_point]['score']) < 0.49:
                continue
            break_points[break_point] = {}
            break_points[break_point]['novel_kmers'] = {}
            break_points[break_point]['singular_kmers'] = {}
            break_points[break_point]['inner_kmers'] = {}
            kmers = track[break_point]['kmers']
            for kmer in kmers:
                count = self.get_kmer_count(kmer, self.index, False)
                if count == 0:
                    break_points[break_point]['novel_kmers'][kmer] = kmers[kmer]
                    n += 1
                # TODO: when to use singular kmers?
                #if count == 1:
                #    break_points[break_point]['singular_kmers'][kmer] = kmers[kmer]
            inner_kmers = track[break_point]['inner_kmers']
            for kmer in inner_kmers:
                count = self.get_kmer_count(kmer, self.index, False)
                if count == 1:
                    break_points[break_point]['inner_kmers'][kmer] = kmers[kmer]
        if n == 0:
            print(red('no novel kmers found for'), white(track_name))
            return None
        path = os.path.join(self.get_current_job_directory(), 'novel_kmers_' + track_name  + '.json') 
        with open(path, 'w') as json_file:
            json.dump({'break_points': break_points}, json_file, sort_keys = True, indent = 4)
        return path

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job to calculate novel kmer overlap
# For the entire set of novel kmers resulting given to as input, finds the structural variations that share each kmer.
# Also outputs for each track some statistics on the level of overlap between it and others
# Plots the overlap statistics and also outputs list of tracks sorted on based on the number of overlapping kmers
# ============================================================================================================================ #
# ============================================================================================================================ #

class NovelKmerOverlapJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = NovelKmerOverlapJob(job_name = 'overlap_', previous_job_name = 'novel_')
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare(self):
        self.minimum_coverage = 5
        with open(os.path.join(self.get_previous_job_directory(), 'merge.json'), 'r') as json_file:
            self.events = json.load(json_file)

    def transform(self, track, track_name):
        novel_kmers = track['novel_kmers']
        overlap = {}
        for event in self.events:
            if event != track_name:
                for kmer in novel_kmers:
                    if novel_kmers[kmer] > self.minimum_coverage:
                        if kmer in self.events[event]['novel_kmers']:
                            if not kmer in overlap:
                                overlap[kmer] = []
                            overlap[kmer].append(event)
        track['overlap'] = overlap
        track['overlap_percentage'] = -1 if len(novel_kmers) == 0 else float(len(overlap)) / float(len(novel_kmers))
        track['overlap_count'] = len(track['novel_kmers']) - len(track['overlap'])
        return track

    def plot(self, tracks):
        self.plot_novel_kmer_overlap_count(tracks)
        self.plot_novel_kmer_overlap_percentage(tracks)

    def plot_novel_kmer_overlap_count(self, tracks):
        path = os.path.join(self.get_current_job_directory(), 'plot_novel_kmer_overlap_count.html')
        x = list(map(lambda x: tracks[x]['overlap_count'], tracks))
        trace = graph_objs.Histogram(
            x = x,
            histnorm = 'count',
            xbins = dict(
                start = 0.0,
                end = 3000.0,
                size = 1.0,
            )
        )
        layout = graph_objs.Layout(title = 'Non-Overlapping Novel kmer Count')
        fig = graph_objs.Figure(data = [trace], layout = layout)
        plotly.plot(fig, filename = path)

    def plot_novel_kmer_overlap_percentage(self, tracks):
        path = os.path.join(self.get_current_job_directory(), 'plot_novel_kmer_overlap_percentage.html')
        x = list(map(lambda x: tracks[x]['overlap_percentage'], tracks))
        trace = graph_objs.Histogram(
            x = x,
            histnorm = 'count',
            xbins = dict(
                start = 0.0,
                end = 1.0,
                size = 0.05
            )
        )
        layout = graph_objs.Layout(title = 'Novel kmer Overlap')
        fig = graph_objs.Figure(data = [trace], layout = layout)
        plotly.plot(fig, filename = path)

    def sort(self, output):
        # can't sort complex dictionary directly, break-up to simpler thing
        tmp = {}
        for track in output:
            tmp[track] = output[track]['overlap_count']
        # sort
        sorted_output = sorted(tmp.items(), key = operator.itemgetter(1))
        # dump
        with open(os.path.join(self.get_current_job_directory(), 'sort.json'), 'w') as json_file:
            json.dump(sorted_output, json_file, sort_keys = True, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce Job to find reads containing kmers from structural variation events that produce too many novel kmers.
# Divides the reads into c.max_threads partitions and counts the kmers inside each partitions
# khmer's results can't be trusted but will help us reduce the exact counting to a target set of kmers.
# Why are we really doing this?
# Only one of the break points for a deletion should happen and that can produce 2*c.ksize but we are seeing way more novel kmers
# than that. We get the reads containing those novel kmers to see what is happenning. For the time being lets only focus on those
# appearing 255 times. Because what the fuck how is a novel kmer appearing 255 times? 
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountKmersExactJob(map_reduce.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Helper to speed up exporting
    # ============================================================================================================================ #

    class ExportHelperJob(map_reduce.Job):

        def find_thread_count(self):
            c = config.Configuration()
            self.num_threads = c.max_threads

        def load_inputs(self):
            for i in range(0, self.num_threads):
                self.batch[i] = {}
            with open(os.path.join(self.get_previous_job_directory(), 'merge.json'), 'r') as json_file:
                paths = json.load(json_file)
                index = 0
                for track in paths:
                    self.batch[index][track] = paths[track]
                    index = index + 1
                    if index == self.num_threads:
                        index = 0

        def transform(self, track, track_name):
            novel_kmers = {}
            with open(track, 'r') as track_file:
                novel_kmers = json.load(track_file)
                for kmer in novel_kmers['novel_kmers']:
                    novel_kmers['novel_kmers'][kmer]['actual_count'] = self.kmers[kmer]
            path = os.path.join(self.get_current_job_directory(), 'exact_counts_' + track_name + '.json')
            with open(path, 'w') as track_file:
                json.dump(novel_kmers, track_file, sort_keys = True, indent = 4)
            return path

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch_export_helper(kmers, **kwargs):
        job = CountKmersExactJob.ExportHelperJob(job_name = 'CountKmersExactJob_', previous_job_name = 'NovelKmerJob_', kmers = kmers, **kwargs)
        job.execute()

    @staticmethod
    def launch(**kwargs):
        job = CountKmersExactJob(job_name = 'CountKmersExactJob_', previous_job_name = 'NovelKmerJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed to specify the set of structural variations (default CHM1_Lumpy.Del.100bp.DEL.bed)
        # --fastq: the genome from which we are getting the kmer counts
        # --thread: number of threads to use
        pass

    def load_inputs(self):
        c = config.Configuration()
        # accumulate all kmers from all tracks into a large dict to improve performance
        # if we keep each event's kmers in a seperate dict, runtime will grow linearly with number of events
        self.kmers = {}
        with open(os.path.join(self.get_previous_job_directory(), 'merge.json'), 'r') as json_file:
            tracks = json.load(json_file)
            for path in tracks:
                with open(tracks[path], 'r') as track_file:
                    track = json.load(track_file)
                    for kmer in track['novel_kmers']:
                        self.kmers[kmer] = 0
        print('counting', len(self.kmers), 'kmers')
        # these two next lines are just to avoid overrding another method
        for index in range(0, self.num_threads):
            self.batch[index] = {}

    def reduce(self):
        c = config.Configuration()
        kmers = self.merge_counts()
        # update counts for each break point
        CountKmersExactJob.launch_export_helper(kmers)
        # export kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as kmers_file:
            json.dump(kmers, kmers_file, sort_keys = True, indent = 4)

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

    def check_cli_arguments(self, args):
        # --bed: the BED file with exon locations
        # --fastq: the sample genome for which we are calculating the coverage
        # --threads: maximum number of threads to use
        # --reference: the reference genome to use for the exon locations (default hg19)
        pass

    def prepare(self):
        self.kmers = {}

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

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
        tracks = bed.read_tracks(c.bed_file)
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

    #def plot(self, _):
        #r = [ self.counts[i] for i in sorted(random.sample(range(len(self.counts)), int(len(self.counts) / 100))) ]
        #data = [graph_objs.Histogram(x = r, xbins = dict(start = 0, end = 500, size = 5))]
        #filename = os.path.join(self.get_current_job_directory(), 'distribution.html')
        #plotly.plot(data, filename = filename, auto_open = False)

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_output_directory(self):
        c = config.Configuration()
        fastq_file_name = c.fastq_file.split('/')[-1][::-1].split('.')[-1][::-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/genotyping/', fastq_file_name, self.job_name[:-1]))

    def get_current_job_directory(self):
        return self.get_output_directory()

# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    if c.job == 'NovelKmersJob':
        NovelKmersJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'CountKmersExactJob':
        CountKmersExactJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'DepthOfCoverageEstimationJob':
        DepthOfCoverageEstimationJob.launch(resume_from_reduce = c.resume_from_reduce)
