import io
import os
import re
import pwd
import sys
import copy
#import json
import math
import time
import random
import argparse
import traceback
import statistics as stats

from kmer import (
    bed,
    sets,
    config,
    commons,
    map_reduce,
    statistics,
    count_server,
)
from kmer.sv import StructuralVariation, Inversion, Deletion
from kmer.kmers import *

from pympler import asizeof

import colorama
import memory_profiler

import rapidjson as json
import plotly.offline as plotly
import plotly.graph_objs as graph_objs

print('done importing')

# ============================================================================================================================ #
# ============================================================================================================================ #
# Once we have extracted the kmers for every breakpoint of a structural variation, we will be interested in knowing
# which of thses kmers are actually novel ones that don't normally occur in the genome.
# REQUIRES a count_server for the reference genome
# ============================================================================================================================ #
# ============================================================================================================================ #

class NovelKmerJob(map_reduce.Job):
    
    # This helper job will read a heavy json file and output each top level object in it to a separate file in parallel
    class NovelKmerLoadInputJob(map_reduce.Job):

        def find_thread_count(self):
            c = config.Configuration()
            self.num_threads = c.max_threads

        def load_inputs(self):
            path = os.path.join(self.get_previous_job_directory(), 'batch_' + str(self.index) + '.json')
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
                for track in self.tracks:
                    path = os.path.join(self.get_current_job_directory(), track)
                    self.tracks[track] = path 

        def transform(self, track, track_name):
            path = os.path.join(self.get_current_job_directory(), track_name)
            print('writing', path)
            with open(path, 'w') as tmp_file:
                json.dump(track, tmp_file, sort_keys = True, indent = 4)
            return path

        def reduce(self):
            # no need to merge stuff here
            pass

        def get_current_job_directory(self):
            # get rid of the final _
            return os.path.abspath(os.path.join(self.get_output_directory(), self.job_name[:-1], 'tracks'))

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = NovelKmerJob(job_name = 'NovelKmerJob_', previous_job_name = 'break_point_')
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
        for index in range(0, self.num_threads):
            print(colorama.Fore.BLUE + 'handling previous batch', index, colorama.Fore.WHITE)
            tracks = self.launch_input_loader(index = index) 
            self.batch[index] = {}
            for track in tracks:
                self.batch[index][track] = tracks[track]

    def launch_input_loader(self, index):
        job = self.NovelKmerLoadInputJob(job_name = 'NovelKmerJob_', previous_job_name = 'break_point_')
        job.index = index
        job.execute()
        return job.tracks

    def transform(self, path, track_name):
        # load each job's input separately to reduce memory footprint
        with open(path, 'r') as json_file:
            track = json.load(json_file)
        remove = {}
        novel_kmers = {}
        for break_point in track:
            # ignore weak break_points
            if float(track[break_point]['score']) < 0.49:
                continue
            kmers = track[break_point]['kmers']
            for kmer in kmers:
                canon = get_canonical_kmer_representation(kmer)
                # we are going to be comparing these againt reads, so it is necessary to consider reverse complement as well
                # khmer returns the same count for a kmer and its reverse complement
                if is_kmer_novel(canon, self.index):
                    if not canon in novel_kmers:
                        novel_kmers[canon] = {
                            'count': [],
                            'break_points': []
                        }
                    # because breakpoint kmers were extracted from reads, thr reverse complement of kmers may not necessarily appear (unless by chance)
                    novel_kmers[canon]['count'] = kmers[kmer]
                    novel_kmers[canon]['break_points'].append(break_point)
        return {'novel_kmers': novel_kmers}

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
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = CountKmersExactJob(job_name = 'CountKmersExactJob_', previous_job_name = 'NovelKmerJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed to specify the set of structural variations
        # --fastq: the genome from which we are getting the kmer counts
        pass

    def load_inputs(self):
        c = config.Configuration()
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        tracks = {}
        with open(path, 'r') as json_file:
            batch = json.load(json_file)
            for track in batch:
                tracks[track] = batch[track]
        # accumulate all kmers from all tracks into a large dict to improve performance
        # if we keep each event's kmers in a seperate dict, runtime will grow linearly with number of events
        self.kmers = {}
        for track in tracks:
            # remember these kmers already come in the canonical representation
            novel_kmers = tracks[track]['novel_kmers']
            for kmer in novel_kmers:
                self.kmers[kmer] = 0
        # we are not gonna need the tracks until the reduce step and it takes a whole lot of memory so delete it
        del tracks
        # these are just dummies
        for index in range(0, self.num_threads):
            self.batch[index] = {} # avoid overrding extra methods from MapReduce

    def reduce(self):
        c = config.Configuration()
        # merge kmer counts from all children
        kmers = self.merge_counts()
        # reindex based on track
        tracks = {}
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            batch = json.load(json_file)
            for track in batch:
                tracks[track] = batch[track]
        for track in tracks:
            novel_kmers = tracks[track]['novel_kmers']
            for novel_kmer in novel_kmers:
                novel_kmers[novel_kmer]['actual_count'] = kmers[novel_kmer]
        # 
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(tracks, json_file, sort_keys = True, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #

class DepthOfCoverageEstimationJob(map_reduce.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = DepthOfCoverageEstimationJob(job_name = 'DepthOfCoverageEstimationJob_', previous_job_name = "", **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # This helper job will extract the set of kmers to be counted for an estimation of depth of coverage
    # As there are many exonic regions to consider and reading them from the reference genome can be time
    # consuming, this job will help save time
    # ============================================================================================================================ #

    class ExtractExonicKmersJob(map_reduce.Job):

        def parse_track_file(self, path):
            tracks = {}
            with open(path) as bed_file:
                line = bed_file.readline()
                while line:
                    tokens = line.split()
                    track = bed.BedTrack(tokens[0], int(tokens[1]), int(tokens[2]))
                    if not track.name in tracks:
                        tracks[track.name] = track
                    line = bed_file.readline()
            return tracks

        # def load_genes(self):
        #     c = config.Configuration()
        #     self.genes = {}
        #     with open(c.genes) as genes_file:
        #         lines = genes_file.readlines()
        #         for line in lines:
        #             gene = line.strip()
        #             print(gene)
        #             self.genes[gene] = True

        def prepare(self):
            self.kmers = {}

        def find_thread_count(self):
            c = config.Configuration()
            self.num_threads = c.max_threads

        def load_inputs(self):
            c = config.Configuration()
            # check if we already have the kmers
            if os.path.isfile(os.path.join(self.get_current_job_directory(), 'kmers.json')):
                with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'r') as kmers_file:
                    print('Exon kmers already extracted, reloading from cache')
                    self.kmers = json.load(kmers_file)
                    # setting num_threads to 0 will bypass all execution
                    self.num_threads = 0
                    return
            # self.load_genes()
            # load exonic regions
            tracks = self.parse_track_file(c.bed_file)
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
            self.kmers[track_name] = {}
            c = config.Configuration()
            seq = bed.extract_sequence(track)
            for kmer in extract_kmers(c.ksize, seq):
                if not kmer in self.kmers[track_name]:
                    self.kmers[track_name][kmer] = 0
                self.kmers[track_name][kmer] += 1
            return None

        def output_batch(self, batch):
            path = os.path.join(self.get_current_job_directory(), 'kmers_' + str(self.index) + '.json')
            with open(path, 'w') as json_file:
                json.dump(self.kmers, json_file, sort_keys = True, indent = 4)

        def reduce(self):
            if self.num_threads == 0:
                return
            # merge all the kmer counts from previous steps
            self.kmers = {}
            for i in range(0, self.num_threads):
                path = os.path.join(self.get_current_job_directory(), 'kmers_' + str(i) + '.json')
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    for track in batch:
                        self.kmers[track] = batch[track]
                        # if not kmer in self.kmers:
                        #     self.kmers[kmer] = 0
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, sort_keys = True, indent = 4)

        def get_output_directory(self):
            c = config.Configuration()
            fastq_file_name = c.fastq_file.split('/')[-1][::-1].split('.')[-1][::-1]
            return os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/genotyping/', fastq_file_name, self.job_name[:-1]))

        def get_current_job_directory(self):
            return self.get_output_directory()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed: the BED file with exon locations
        # --fastq: the sample genome for which we are calculating the coverage
        # --threads: maximum number of threads to use
        # --reference: the reference genome to use for the exon locations (default hg19)
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        c = config.Configuration()
        print('running helper job to get kmers')
        job = self.ExtractExonicKmersJob(job_name = 'DepthOfCoverageEstimationJob_', previous_job_name = '_', resume_from_reduce = False)
        job.execute()
        # all kmers begin with a count of zero
        tracks = job.kmers
        kmers = {}
        for track in tracks:
            for kmer in tracks[track]:
                kmers[kmer] = 0
        self.kmers = {}
        # sample from these kmers
        print('got', len(kmers), 'kmers, sampling 100000 randomly')
        sample = random.sample(kmers.items(), 100000)
        for key, value in sample:
            self.kmers[key] = value
        print('kmers available, counting for coverage')
        # dummy, to avoid overrding distribute_workload
        for index in range(0, self.num_threads):
            self.batch[index] = {}

    def reduce(self):
        c = config.Configuration()
        # merge kmer counts from all children
        kmers = self.merge_counts()
        # calculate mean and std
        self.counts = list(map(lambda x: kmers[x], list(kmers.keys())))
        self.mean = stats.mean(self.counts)
        # filter outliers
        self.counts = list(filter(lambda x: x < 3 * self.mean, self.counts))
        self.mean = stats.mean(self.counts)
        # filter anything appearing more than twice the medium, 4x coverage or more, repeatimg kmer
        self.counts = list(filter(lambda x: x < 2 * self.mean, self.counts))
        self.mean = stats.mean(self.counts)
        self.std = stats.stdev(self.counts)
        print('mean:', self.mean)
        print('std:', self.std)
        kmers = {'kmers': kmers}
        kmers['coverage'] = {'mean': self.mean, 'std': self.std}
        # output merged kmer counts
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4)
        self.plot(self.counts)

    def plot(self, _):
        r = [ self.counts[i] for i in sorted(random.sample(range(len(self.counts)), int(len(self.counts) / 100))) ]
        data = [graph_objs.Histogram(x = r, xbins = dict(start = 0, end = 500, size = 5))]
        filename = os.path.join(self.get_current_job_directory(), 'distribution.html')
        plotly.plot(data, filename = filename, auto_open = False)

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
    #CountKmersExactJob.launch()
    #NovelKmerJob.launch()
    DepthOfCoverageEstimationJob.launch(resume_from_reduce = False)
