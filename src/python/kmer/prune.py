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
    count_server,
)
from kmer.threading import (
    exact_counter
)
from kmer.sv import StructuralVariation, Inversion, Deletion

import colorama
import memory_profiler
import plotly.offline as plotly
import plotly.graph_objs as graph_objs

from pympler import asizeof

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

def get_canonical_kmer_representation(kmer):
    reverse_complement = bed.reverse_complement_sequence(kmer)
    return kmer if kmer < reverse_complement else reverse_complement

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
            self.prefix = 'tmp_' + str(self.index) + '_'
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
                    path = os.path.join(self.get_current_job_directory(), self.prefix + track)
                    self.tracks[track] = path 

        def transform(self, track, track_name):
            path = os.path.join(self.get_current_job_directory(), self.prefix + track_name)
            print('writing', path)
            with open(path, 'w') as tmp_file:
                json.dump(track, tmp_file, sort_keys = True, indent = 4)
            return path

        def reduce(self):
            # no need to merge stuff here
            pass

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
    #TODO: test this
    def load_inputs(self):
        for index in range(0, self.num_threads):
            print(colorama.Fore.BLUE + 'handling previous batch', index, colorama.Fore.WHITE)
            tracks = self.launch_input_loader(index = index) 
            self.batch[index] = {}
            for track in tracks:
                self.batch[index][track] = tracks[track]
            #path = os.path.join(self.get_previous_job_directory(), 'batch_' + str(index) + '.json')
            #with open(path, 'r') as json_file:
            #    batch = json.load(json_file)
            #    for track in batch:
            #        # no need to write these files again if we have created them previously
            #        path = os.path.join(self.get_current_job_directory(), 'tmp_' + str(index) + '_' + track)
            #        with open(path, 'w') as tmp_file:
            #            json.dump(batch[track], tmp_file, sort_keys = True, indent = 4)
            #        self.batch[index][track] = path

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
            json.dump(sorted_output, json_file, sort_keys = True, indent = 4, separators = (',', ': '))

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
    def launch(**kwargs):
        job = CountKmersExactJob(job_name = 'CountKmersExactJob_', previous_job_name = 'NovelKmerJob_', **kwargs)
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
        ahead = self.fastq_file.readline()
        n = 0
        m = 0
        t = time.time()
        while ahead:
            if state == HEADER_LINE:
                if line[0] == '@' and ahead[0] != '@':
                    if self.fastq_file.tell() >= (self.index + 1) * self.fastq_file_chunk_size:
                        print(self.index, 'reached segment boundary')
                        break
                    state = SEQUENCE_LINE
                    name = line[:-1] # ignore the EOL character
            elif state == SEQUENCE_LINE:
                state = THIRD_LINE
                seq = line[:-1] # ignore the EOL character
                n += 1
                if n == 100000:
                    n = 0
                    m += 1
                    c = self.fastq_file.tell() - self.index * self.fastq_file_chunk_size
                    s = time.time()
                    p = c / float(self.fastq_file_chunk_size)
                    e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                    print(self.index, 'progress:', p, 'took: ', s - t, 'ETA: ', e)
                    #print(self.index, 'm =', m)
                yield seq, name
            elif state == THIRD_LINE:
                state = QUALITY_LINE
            elif state == QUALITY_LINE:
                state = HEADER_LINE
            line = ahead
            ahead = self.fastq_file.readline()
        print(self.index, ' end of input')

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --fastq option to specify the file to read from
        pass

    def prepare(self):
        self.minimum_coverage = 5

    def load_inputs(self):
        c = config.Configuration()
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        tracks = {}
        with open(path, 'r') as json_file:
            batch = json.load(json_file)
            for track in batch:
                tracks[track] = batch[track]
        #print('tracks:', asizeof.asizeof(tracks))
        self.kmers = {}
        for track in tracks:
            # remember these kmers already come in the canonical representation
            novel_kmers = tracks[track]['novel_kmers']
            for kmer in novel_kmers:
                self.kmers[kmer] = 0
        del tracks
        #print('kmers:', asizeof.asizeof(self.kmers))
        for index in range(0, self.num_threads):
            self.batch[index] = {} # avoid overrding extra methods from MapReduce

    def run_batch(self, batch):
        c = config.Configuration()
        self.fastq_file = open(c.fastq_file, 'r')
        self.fastq_file_chunk_size = math.ceil(os.path.getsize(self.fastq_file.name) / float(self.num_threads))
        self.fastq_file.seek(self.index * self.fastq_file_chunk_size, 0)
        # this forked process will exit at the end of the following function call
        self.transform()
        self.output_batch(self.kmers)

    def transform(self):
        c = config.Configuration()
        # it is gonna take too long to search inside each track when parsing the FASTQ file so lets
        # re-index based on kmer rather than track (1 table lookup instead of len(tracks))
        #output = {}
        #for track in self.tracks:
        #    # remember these kmers already come in the canonical representation
        #    novel_kmers = self.tracks[track]['novel_kmers']
        #    for kmer in novel_kmers:
        #        output[kmer] = 0
        #
        for read, name in self.parse_fastq():
            kmers = get_all_kmers(read, c.ksize)
            for kmer in kmers:
                # novel kmers for each track are already stored in canonical representation
                canon = get_canonical_kmer_representation(kmer)
                if canon in self.kmers: 
                    self.kmers[canon] += 1

    def reduce(self):
        c = config.Configuration()
        # merge kmer counts from all children
        kmers = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            if not os.path.isfile(path):
                print(colorama.Fore.RED + 'couldn\'t find batch', i, ' results will be suspicious')
                continue
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    if not kmer in kmers:
                        kmers[kmer] = 0
                    kmers[kmer] += batch[kmer]
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
            json.dump(self.tracks, json_file, sort_keys = True, indent = 4, separators = (',', ': '))

class CountKmersExactMultiThread():

    @staticmethod
    def launch():
        exact_counter.execute()

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
        job = CountBedKmersExactJob(job_name = 'CountBedKmersExactJob_', previous_job_name = 'KmerNormalDistributionFittingJob_', resume_from_reduce = True)
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
        tmp = {}
        with open(path, 'r') as json_file:
            tmp = json.load(json_file)
        kmers = random.sample(tmp.items(), int(len(tmp) / 10))
        for pair in kmers:
            self.kmers[pair[0]] = pair[1]
        print(len(self.kmers))
        print('done loading inputs')

    def transform(self):
        print('transforming', self.index)
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
        return output

    def reduce(self):
        c = config.Configuration()
        # 
        output = {}
        for i in range(0, self.num_threads):
            print('batch', i)
            try:
                with open(os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json'), 'r') as json_file:
                    batch = json.load(json_file)
                    for kmer in batch:
                        reverse = bed.reverse_complement_sequence(kmer)
                        if kmer in output:
                            output[kmer] += batch[kmer]
                        elif reverse in output:
                            output[reverse] += batch[kmer] 
                        else:
                            output[kmer] = batch[kmer]
            except:
                pass
        #
        self.counts = list(map(lambda x: output[x], list(output.keys())))
        self.mean = stats.mean(self.counts)
        # filter outliers
        self.counts = list(filter(lambda x: x < 3 * self.mean, self.counts))
        self.mean = stats.mean(self.counts)
        self.std = stats.stdev(self.counts)
        print('mean:', self.mean)
        print('std:', self.std)
        # 
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        self.plot()

    def plot(self):
        n = statistics.NormalDistribution(mean = self.mean, std = self.std)
        r = [ self.counts[i] for i in sorted(random.sample(range(len(self.counts)), int(len(self.counts) / 100))) ]
        #y = list(map(lambda x: n.pmf(x), r))
        #trace = graph_objs.Scatter(x = r, y = y, mode = 'markers')
        #data = [trace]
        data = [graph_objs.Histogram(x = r, xbins = dict(start = 0, end = 500, size = 5))]
        path = os.path.join(self.get_current_job_directory(), 'distribution')
        plotly.plot(data, filename = path + '.html', auto_open = False)

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
    #NovelKmerJob.launch()
    CountKmersExactJob.launch(resume_from_reduce = True)
