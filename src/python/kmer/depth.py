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
    counter,
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
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExonDepthOfCoverageCountingJob(counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ExonDepthOfCoverageCountingJob(job_name = 'ExonDepthOfCoverageCountingJob_', previous_job_name = "", category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        # check if we already have the kmers
        if os.path.isfile(os.path.join(self.get_current_job_directory(), 'kmers.json')):
            with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'r') as kmers_file:
                print('Exon kmers already extracted, reloading from cache')
                self.kmers = json.load(kmers_file)
                # setting num_threads to 0 will bypass all execution and jump to reduce
                self.num_threads = 0
                exit()
        # load exonic regions
        tracks = bed.read_tracks(c.exons)
        for track in tracks:
            seq = bed.extract_sequence(tracks[track])
            for kmer in extract_canonical_kmers(c.ksize, seq):
                self.kmers[kmer] = {'count': 0}
        # split exonic regions into several batches
        for i in range(0, self.num_threads):
            self.batch[i] = {}

    def reduce(self):
        self.merge_counts()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExonDepthOfCoverageEstimationJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = DepthOfCoverageEstimationJob(job_name = 'DepthOfCoverageEstimationJob_', previous_job_name = 'DepthOfCoverageCountingJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare(self):
        self.kmers = {}

    def load_inputs(self):
        c = config.Configuration()
        #self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        self.counts_provider = counttable.DictionaryCountsProvider(json.load(open(os.path.join(self.get_previous_job_directory(), 'kmers.json'), 'r')))
        # check if we already have the kmers
        #if os.path.isfile(os.path.join(self.get_current_job_directory(), 'kmers.json')):
        #    with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'r') as kmers_file:
        #        print('Exon kmers already extracted, reloading from cache')
        #        self.kmers = json.load(kmers_file)
        #        # setting num_threads to 0 will bypass all execution and jump to reduce
        #        self.num_threads = 0
        #        return
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
        for kmer in extract_canonical_kmers(c.ksize, seq):
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
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueKmersDepthOfCoverageEstimationJob(map_reduce.BaseGenotypingJob, counter.SimulationExactCountingJob):# if config.Configuration().simulation else counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = UniqueKmersDepthOfCoverageEstimationJob(job_name = 'UniqueKmersDepthOfCoverageEstimationJob_', previous_job_name = 'UniqueKmersDepthOfCoverageCountingJob_', category = 'programming', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        print(self.get_current_job_directory())
        c = config.Configuration()
        #self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[1])
        self.kmers = {}
        n = 0
        self.acc = {}
        for i in range(1, 2):
            self.acc[i] = 0
        for kmer, count in self.reference_counts_provider.stream_kmers():
            if count in self.acc:
                if self.acc[count] == 1000000:
                    self.acc.pop(count, None)
                    if len(self.acc) == 0:
                        break
                    continue
                self.kmers[kmer] = {'ref': count, 'count': 0}
                self.kmers[reverse_complement(kmer)] = {'ref': count, 'count': 0}
                self.acc[count] += 1
                if self.acc[count] % 1000 == 0:
                    print(self.acc[count], 'kmers with a count of', count, 'so far')
        print('Counting', green(len(self.kmers)), 'unique kmers')
        self.round_robin()

    def reduce(self):
        self.kmers = self.merge_counts()
        self.plot_reference_distribution(self.kmers)
        self.counts = list(map(lambda kmer: self.kmers[kmer]['count'], self.kmers))
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

    def plot_reference_distribution(self, kmers):
        r = []
        for i in range(1, 2):
            k = list(filter(lambda kmer: kmers[kmer]['ref'] == i, kmers))
            counts = list(map(lambda kmer: kmers[kmer]['count'], k))
            mean = numpy.mean(counts)
            counts = list(filter(lambda x: x < 3 * mean, counts))
            mean = numpy.mean(counts)
            counts = list(filter(lambda x: x < 3 * mean, counts))
            r.append(counts)
        visualizer.bar(ys = [list(map(lambda x: statistics.mean(x), r)), list(map(lambda x: statistics.variance(x), r))], x = list(range(1, 2)), name = 'bar plot mean kmer coverage by reference frequency', path = self.get_current_job_directory(), x_label = 'reference frequency', y_label = 'mean coverage')

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_current_job_directory(self):
        return map_reduce.Job.get_current_job_directory(self)

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
