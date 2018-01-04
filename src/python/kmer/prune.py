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
import statistics as stats

from kmer import (
    bed,
    sets,
    config,
    commons,
    counttable,
    map_reduce,
    statistics,
)
from kmer.sv import StructuralVariation, Inversion, Deletion

print('importing khmer')
#import khmer
import colorama
#import pybedtools
import plotly.offline as plotly
import plotly.graph_objs as graph_objs

print('done importing')
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
# MapReduce to find reads containing kmers from structural variation events that produce too many novel kmers.
# khmer's results simply can't be trusted but will help us reduce the exact counting to a target set of kmers.
# Why are we really doing this?
# Only one of the break points for a deletion should happen and that can produce 2*c.ksize but we are seeing way more novel kmers
# than that. We get the reads containing those novel kmers to see what is happenning. For the time being lets only focus on those
# appearing 255 times. Because what the fuck how is a novel kmer appearing 255 times? 
# ============================================================================================================================ #
# ============================================================================================================================ #

def get_all_kmers(read, k):
    kmers = []
    for i in range(0, len(read) - k + 1):
        kmer = read[i : i + k]
        kmers.append(kmer)
    return kmers

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
        m = 0
        t = time.time()
        while line:
            if state == HEADER_LINE:
                if line[0] == '@':
                    if self.fastq_file.tell() >= (self.index + 1) * self.fastq_file_chunk_size:
                        print(self.index, 'reached segment boundary')
                        break
                    state = SEQUENCE_LINE
                    name = line[:-1] # ignore the EOL character
                line = self.fastq_file.readline()
                continue
            if state == SEQUENCE_LINE:
                #if m > 47:
                #    print(self.index, 'SEQ', n, m)
                state = THIRD_LINE
                seq = line[:-1] # ignore the EOL character
                n += 1
                if n == 1000:
                    n = 0
                    m += 1
                    c = self.fastq_file.tell() - self.index * self.fastq_file_chunk_size
                    s = time.time()
                    p = c / float(self.fastq_file_chunk_size)
                    e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                    #print(self.index, 'progress:', p, 'took: ', s - t, 'ETA: ', e)
                    #print(self.index, 'm =', m)
                if m == 50:
                    break
                #if m > 47:
                #    print(self.index, 'yil', n, m)
                yield seq, name
                #if m > 47:
                #    print(self.index, 'red', n, m)
                line = self.fastq_file.readline()
                #if m > 47:
                #    print(self.index, 'cnt', n, m)
                continue
            if state == THIRD_LINE:
                state = QUALITY_LINE
                line = self.fastq_file.readline()
                continue
            if state == QUALITY_LINE:
                state = HEADER_LINE
                line = self.fastq_file.readline()
                continue
        print(self.index, ' end of input')

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
        # TODO: why not load the 'merge.json' and search that one instead?
        for index in range(0, self.num_threads):
            self.batch[index] = {} # avoid overrding extra methods from MapReduce
            path = os.path.join(self.get_previous_job_directory(), 'batch_' + str(index) + '.json')
            with open(path, 'r') as json_file:
                batch = json.load(json_file)
                # find the tracks we are interested, each one might be in a different batch
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
        # this forked process will exit at the end of the following function call 
        self.output_batch(self.transform())

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
            kmers = get_all_kmers(read, c.ksize)
            for kmer in kmers:
                # consider reverse complement as well
                reverse_complement = bed.reverse_complement_sequence(kmer)
                for track in self.tracks:
                    novel_kmers = self.tracks[track]['novel_kmers']
                    if kmer in novel_kmers or reverse_complement in novel_kmers:
                        normal = False
                        reverse = False
                        if kmer in novel_kmers:
                            novel_kmer = kmer
                            normal = True
                        if reverse_complement in novel_kmers:
                            novel_kmer = reverse_complement
                            reverse = True
                        if normal and reverse:
                            print('BOTH APPEAR')
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
# ============================================================================================================================ #
# Like CountKmersExactJob but it reads the kmer it needs to count from the output of a previous ExtractBedKmersJob
# ============================================================================================================================ #
# ============================================================================================================================ #

class CountBedKmersExactJob(CountKmersExactJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = CountBedKmersExactJob(job_name = 'CountBedKmersExactJob_', previous_job_name = 'KmerNormalDistributionFittingJob_')
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --bed: name of the bed file the kmers came from, needs to match the value used for KmerNormalDistributionFittingJob
        # --fastq: path to the fastq file to count these kmers in
        pass

    def prepare(self):
        self.kmers = {}

    def load_inputs(self):
        print('loading inputs...')
        c = config.Configuration()
        # 
        for index in range(0, self.num_threads):
            self.batch[index] = {} # avoid overrding extra methods from MapReduce
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            self.kmers = json.load(json_file)
        print(len(self.kmers))
        print('done loading inputs')

    def transform(self):
        #print('transforming', self.index)
        c = config.Configuration()
        # this one is more complicated because most results are to be commited to a global data store
        # create a copy of data for this instance
        output = {}
        # 
        for read, name in self.parse_fastq():
            kmers = get_all_kmers(read, c.ksize)
            for kmer in kmers:
                reverse_complement = bed.reverse_complement_sequence(kmer)
                if kmer in self.kmers:
                    if not kmer in output:
                        output[kmer] = 0
                    output[kmer] += 1
                if reverse_complement in self.kmers:
                    if not kmer in output:
                        output[kmer] = 0
                    output[kmer] += 1
        print('transformation done ', self.index)
        with open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.out'), 'w') as output_file:
            output_file.write('tranform')
        return output

    def reduce(self):
        c = config.Configuration()
        # 
        output = {}
        for i in range(0, self.num_threads):
            print('batch', i)
            with open(os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json'), 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    if not kmer in output:
                        output[kmer] = 0
                    output[kmer] += batch[kmer]
        #
        self.counts = list(map(lambda x: output[x], list(output.keys())))
        self.median = stats.median(self.counts)
        self.std = stats.stdev(self.counts)
        print('median:', self.median)
        print('std:', self.std)
        # 
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        self.plot()

    def plot(self):
        n = statistics.NormalDistribution(mean = self.median, std = self.std)
        r = [ self.counts[i] for i in sorted(random.sample(range(len(self.counts)), int(len(self.counts) / 50))) ]
        y = list(map(lambda x: n.pmf(x), r))
        trace = graph_objs.Scatter(x = r, y = y, mode = 'markers')
        data = [trace]
        path = os.path.join(self.get_current_job_directory(), 'distribution')
        #layout = graph_objs.Layout(title = 'Exon kmer Coutn Distribtuion')
        #fig = graph_objs.Figure(data = data, layout = layout)
        plotly.plot(data, filename = path + '.html', auto_open = False, image = 'png', image_filename = path + '.png')
        #plotly.image.save_as(fig, filename = path)

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractBedKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # job-specific stuff
    # ============================================================================================================================ #

    def parse_track_file(self, path):
        tracks = {}
        with open(path) as bed_file:
            line = bed_file.readline()
            while line:
                tokens = line.split()
                track = bed.BedTrack(tokens[-3], int(tokens[-2]), int(tokens[-1]))
                if not track.name in tracks:
                    tracks[track.name] = track
                line = bed_file.readline()
        return tracks

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = KmerNormalDistributionFittingJob(job_name = 'KmerNormalDistributionFittingJob_', previous_job_name = "", **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare(self):
        self.kmers = {}

    def check_cli_arguments(self, args):
        # --reference: the reference genome to use
        # --bed: the BED file to use
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        c = config.Configuration()
        # 
        tracks = self.parse_track_file(c.bed_file)
        for i in range(0, self.num_threads):
            self.batch[i] = {} # avoid overrding extra methods from MapReduce
        #
        index = 0
        for track in tracks:
            self.batch[index][track] = tracks[track]
            index += 1
            if index == self.num_threads:
                index = 0

    def transform(self, track, track_name):
        c = config.Configuration()
        seq = bed.extract_sequence(track)
        for kmer in get_all_kmers(seq, c.ksize):
            if not kmer in self.kmers:
                self.kmers[kmer] = 0
            self.kmers[kmer] += 1
        # we are not chaning the input track, storing results in a object property instead
        return track

    def output_batch(self, batch):
        with open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        exit()

    def reduce(self):
        c = config.Configuration()
        kmers = {}
        # merge all the kmer counts from previous steps
        for i in range(0, self.num_threads):
            print('batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    for kmer in batch:
                        if not kmer in kmers:
                            kmers[kmer] = 0
                        kmers[kmer] += batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4, separators = (',', ': '))

# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    print('main')
    config.init()
    # 
    CountBedKmersExactJob.launch()
