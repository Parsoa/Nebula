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
import traceback

from multiprocessing import Process

from kmer import (
    bed,
    sets,
    config,
    commons,
    counttable,
    map_reduce,
    break_point,
    count_server,
)
from kmer.sv import StructuralVariation, Inversion, Deletion

import khmer
import colorama
import pybedtools
import plotly.offline as plotly
import plotly.graph_objs as graph_objs

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

def get_novel_kmers(kmers, index):
    novel_kmers = {}
    for kmer in kmers:
        if not kmer in novel_kmers:
            count = count_server.get_kmer_count(kmer, index, True)
            if count != 0:
                # this is a novel kmer
                novel_kmers[kmer] = True
    return novel_kmers

def has_novel_kmers(kmers, index):
    # checks if this candidate has a kmer that has not occured in the reference genome
    for kmer in kmers:
        if is_kmer_novel(kmer, index):
            return True
    return False

def is_kmer_novel(kmer, index):
    count = count_server.get_kmer_count(kmer, index, True)
    return count == 0

def has_unique_novel_kmers(track, candidate, kmers, index):
    # checks if this candidate has a novel kmer that hasn't occurred in any other candidate
    for kmer in kmers:
        if is_kmer_novel(kmer, index):
            found = False
            for break_point in track:
                # skip the wicked candidate count key
                if break_point.find('candidates') != -1:
                    continue
                if candidate != break_point:
                    # this kmer appears in at least one other break point so no need to go further
                    if kmer in track[break_point]['kmers']:
                        found = True
                        break
            # we didn't find this novel kmer anywhere so it should be unique, no need to check others
            if not found:
                return True
    # we haven't returned yet so we didn't find any uniquely novel kmers
    return False

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job to extracting high-coverage novel kmers
# ============================================================================================================================ #
# ============================================================================================================================ #

class NovelKmerJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = NovelKmerJob(job_name = 'novel_', previous_job_name = 'break_point_')
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        remove = {}
        novel_kmers = {}
        for candidate in track:
            # remove all candidates eventually
            remove[candidate] = True
            # skip the json key holding the number of candidates
            if candidate.find('candidates') != -1:
                continue
            kmers = track[candidate]['kmers']
            for kmer in kmers:
                if is_kmer_novel(kmer, self.index):
                    if not kmer in novel_kmers:
                        novel_kmers[kmer] = {
                            'count': [],
                            'break_points': []
                        }
                    novel_kmers[kmer]['count'] = kmers[kmer]
                    novel_kmers[kmer]['break_points'].append(candidate)
        # remove every candidate, make oupput more concise
        for candidate in remove:
            track.pop(candidate, None)
        # keep only those with high enough coverage
        track['novel_kmers'] = novel_kmers
        return track

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
            json.dump(sorted_output, json_file, sort_keys = True, indent = 4, separators = (',', ': '))

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce to find reads containing kmers from deletion that produce too many novel kmers.
# Why are we really doing this?
# Only one of the break points for a deletion should happen and that can produce 2*c.ksize but we are seeing way more novel kmers
# than that. We get the reads containing those novel kmers to see what is happenning. For the time being lets only focus on those
# appearing 255 times. Because what the fuck how is a novel kmer appearing 255 times? 
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountKmersExactJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = CountKmersExactJob(job_name = 'exact_', previous_job_name = 'novel_')
        job.execute()

    # ============================================================================================================================ #
    # job-specific stuff
    # ============================================================================================================================ #

    def parse_fastq(self):
        name = None
        HEADER_LINE = 0
        SEQUENCE_LINE = 1
        THIRD_LINE = 2
        QUALITY_LINE = 3
        state = HEADER_LINE
        # need to skip invalid lines
        line = self.fastq_file.readline()
        n = 0
        t = time.time()
        while line:
            if state == HEADER_LINE:
                if line[0] == '@':
                    if self.fastq_file.tell() >= (self.index + 1) * self.fastq_file_chunk_size:
                        break
                    state = SEQUENCE_LINE
                    name = line[:-1] # ignore the EOL character
                line = self.fastq_file.readline()
                continue
            if state == SEQUENCE_LINE:
                state = THIRD_LINE
                seq = line[:-1] # ignore the EOL character
                n += 1
                if n == 5000:
                    n = 0
                    c = self.fastq_file.tell() - self.index * self.fastq_file_chunk_size
                    s = time.time()
                    p = c / float(self.fastq_file_chunk_size)
                    e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                    print(self.index, 'progress:', p, 'took: ', s - t, 'ETA: ', e)
                yield seq, name
                line = self.fastq_file.readline()
                continue
            if state == THIRD_LINE:
                state = QUALITY_LINE
                line = self.fastq_file.readline()
                continue
            if state == QUALITY_LINE:
                state = HEADER_LINE
                line = self.fastq_file.readline()
                continue

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # requires only --fastq option to specify the file to read from
        pass

    def prepare(self):
        self.event_names = ["chr5_78277739_78278045", "chr6_51256300_51256602",\
            "chr18_75997296_75997615", "chr11_102471858_102472162", "chr2_57342743_57343056"]
        self.minimum_coverage = 5
        self.tracks = {}

    def load_inputs(self):
        c = config.Configuration()
        # 
        for index in range(0, self.num_threads):
            self.batch[index] = {} # avoid overrding extra methods from MapReduce
            path = os.path.join(self.get_previous_job_directory(), 'batch_' + str(index) + '.json')
            with open(path, 'r') as json_file:
                batch = json.load(json_file)
                # find the tracks we are interested
                for track in batch:
                    for event_name in self.event_names:
                        if event_name in track:
                            self.tracks[event_name] = batch[track]
                            print('found', track)

    def run_batch(self, batch):
        c = config.Configuration()
        self.fastq_file = open(c.fastq_file, 'r')
        self.fastq_file_chunk_size = math.ceil(os.path.getsize(self.fastq_file.name) / float(self.num_threads))
        self.fastq_file.seek(self.index * self.fastq_file_chunk_size, 0)
        # 
        self.output_batch(self.transform())
        # 
        print(colorama.Fore.GREEN + 'process ', self.index, ' done')

    def transform(self):
        c = config.Configuration()
        # this one is more complicated because most results are to be commited to a global data store
        # create a copy of data for this instance
        output = {}
        for track in self.tracks:
            output[track] = {
                'novel_kmers': {}
            }
        # 
        for read, name in self.parse_fastq():
            for track in self.tracks:
                novel_kmers = self.tracks[track]['novel_kmers']
                for novel_kmer in novel_kmers:
                    # only look at those which appear more than a threshold
                    if novel_kmers[novel_kmer]['count'] >= self.minimum_coverage:
                        # consider reverse complement as well
                        reverse_complement = bed.reverse_complement_sequence(novel_kmer)
                        # read supports this novel kmer
                        if read.find(novel_kmer) != -1 or read.find(reverse_complement) != -1:
                            # add the kmer to output only now, should reduce memory footprint
                            if not novel_kmer in output[track]['novel_kmers']:
                                output[track]['novel_kmers'][novel_kmer] = {
                                    'count': novel_kmers[novel_kmer]['count'],
                                    'break_points': novel_kmers[novel_kmer]['break_points'],
                                    'actual_count': 0,
                                }
                            # update the exact count of the kmer
                            output[track]['novel_kmers'][novel_kmer]['actual_count'] += 1
        return output

    def reduce(self):
        c = config.Configuration()
        # 
        output = {}
        for track in self.tracks:
            output[track] = {
                'novel_kmers': {}
            }
        # 
        for i in range(0, self.num_threads):
            with open(os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json'), 'r') as json_file:
                batch = json.load(json_file)
                for track in batch:
                    novel_kmers = batch[track]['novel_kmers']
                    for novel_kmer in novel_kmers:
                        if not novel_kmer in output[track]['novel_kmers']:
                            output[track]['novel_kmers'][novel_kmer] = {
                                'actual_count': 0
                            }
                        output[track]['novel_kmers'][novel_kmer]['actual_count'] += novel_kmers[novel_kmer]['actual_count']
                        output[track]['novel_kmers'][novel_kmer]['break_points'] = novel_kmers[novel_kmer]['break_points']
                        output[track]['novel_kmers'][novel_kmer]['count'] = novel_kmers[novel_kmer]['count']
        # 
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4, separators = (',', ': '))

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    # 
    CountKmersExactJob.launch()
