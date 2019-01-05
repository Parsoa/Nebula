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
import subprocess

from kmer import (
    bed,
    config,
    counter,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.debug import *
from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

import numpy

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueKmersDepthOfCoverageEstimationJob(map_reduce.BaseGenotypingJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'UniqueKmersDepthOfCoverageEstimationJob'
    _category = 'programming'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = UniqueKmersDepthOfCoverageEstimationJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        print(self.get_current_job_directory())
        c = config.Configuration()
        self.load_reference_counts_provider() 
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
                self.kmers[canonicalize(kmer)] = {'ref': count, 'count': 0}
                self.acc[count] += 1
                if self.acc[count] % 1000 == 0:
                    print(self.acc[count], 'kmers with a count of', count, 'so far')
        print('Counting', green(len(self.kmers)), 'unique kmers')
        if c.accelerate:
            with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4)
        self.round_robin()

    def merge_accelerator_counts(self):
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
            n = 0
            with open (path, 'r') as json_file:
                line = json_file.readline()
                while line:
                    try:
                        kmer = line[:line.find(':')]
                        i = line.find(':')
                        count = int(line[i + 1:])
                        canon = canonicalize(kmer)
                        self.kmers[canon]['count'] += count / 2
                        line = json_file.readline()
                        n += 1
                    except Exception as e:
                        print(n, i, line)
                        print(e)
                        traceback.print_exc()
                        debug_breakpoint()
                        line = json_file.readline()
                        n += 1

    def reduce(self):
        #self.kmers = self.merge_counts()
        self.merge_accelerator_counts()
        self.counts = list(map(lambda kmer: self.kmers[kmer]['count'], self.kmers))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter outliers
        self.counts = list(filter(lambda x: x < 3 * self.mean, self.counts))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter anything appearing more than twice the medium, 4x coverage or more, repeatimg kmer
        self.counts = list(filter(lambda x: x < 2 * self.mean, self.counts))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        #
        self.plot_reference_distribution([ self.counts[i] for i in sorted(random.sample(xrange(len(self.counts)), 10000)) ])
        with open(os.path.join(self.get_current_job_directory(), 'stats_' + str(c.ksize) + '.json'), 'w') as json_file:
            json.dump({ 'mean': self.mean, 'std': self.std }, json_file, sort_keys = True, indent = 4)

    def plot_reference_distribution(self, counts):
        visualizer.histogram(counts, name = 'kmer count distribution', path = self.get_current_job_directory(), x_label = 'count', y_label = 'number of kmers')

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_current_job_directory(self):
        c = config.Configuration()
        if c.simulation:
            return map_reduce.Job.get_current_job_directory(self)
        else:
            bed_file_name = c.bed_file.split('/')[-1]
            return os.path.abspath(os.path.join(self.get_output_directory(), self._name))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class GappedKmersDepthOfCoverageEstimationJob(map_reduce.BaseGenotypingJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'GappedKmersDepthOfCoverageEstimationJob'
    _category = 'programming'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = GappedKmersDepthOfCoverageEstimationJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #
    # ============================================================================================================================ #

    class ExtractExonGappedKmersHelper(map_reduce.BaseGenotypingJob):

        _name = 'ExtractExonGappedKmersHelper'
        _category = 'programming'
        _previous_job = None

        @staticmethod
        def launch(**kwargs):
            job = ExtractExonGappedKmersHelper(**kwargs)
            job.execute()
        
        def load_inputs(self):
            c = config.Configuration()
            exons = bed.load_tracks_from_file_as_dict(c.exons)
            extract_whole_genome()
            self.round_robin(exons)
            self.kmers = {}

        def transform(self, track, track_name):
            chromosome = extract_chromosome(self.chrom.lower())
            sequence = chromosome[track.begin: track.end]
            for i in range(0, len(sequence) - 50):
                kmer = canonicalize(sequence(i, 37))
                self.kmers[kmer] = 0
            return None

        def output_batch(self, batch):
            json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
            json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
            json_file.close()
            exit()

        def reduce(self):



    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        print(self.get_current_job_directory())
        c = config.Configuration()
        self.load_reference_counts_provider() 
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
                self.kmers[canonicalize(kmer)] = {'ref': count, 'count': 0}
                self.acc[count] += 1
                if self.acc[count] % 1000 == 0:
                    print(self.acc[count], 'kmers with a count of', count, 'so far')
        print('Counting', green(len(self.kmers)), 'unique kmers')
        if c.accelerate:
            with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4)
        self.round_robin()

    def merge_accelerator_counts(self):
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
            n = 0
            with open (path, 'r') as json_file:
                line = json_file.readline()
                while line:
                    try:
                        kmer = line[:line.find(':')]
                        i = line.find(':')
                        count = int(line[i + 1:])
                        canon = canonicalize(kmer)
                        self.kmers[canon]['count'] += count / 2
                        line = json_file.readline()
                        n += 1
                    except Exception as e:
                        print(n, i, line)
                        print(e)
                        traceback.print_exc()
                        debug_breakpoint()
                        line = json_file.readline()
                        n += 1

    def reduce(self):
        #self.kmers = self.merge_counts()
        self.merge_accelerator_counts()
        self.counts = list(map(lambda kmer: self.kmers[kmer]['count'], self.kmers))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter outliers
        self.counts = list(filter(lambda x: x < 3 * self.mean, self.counts))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter anything appearing more than twice the medium, 4x coverage or more, repeatimg kmer
        self.counts = list(filter(lambda x: x < 2 * self.mean, self.counts))
        self.mean = numpy.mean(self.counts)
        self.std = numpy.std(self.counts)
        print(len(self.counts))
        print('mean:', self.mean)
        print('std:', self.std)
        #
        self.plot_reference_distribution([ self.counts[i] for i in sorted(random.sample(xrange(len(self.counts)), 10000)) ])
        with open(os.path.join(self.get_current_job_directory(), 'stats_' + str(c.ksize) + '.json'), 'w') as json_file:
            json.dump({ 'mean': self.mean, 'std': self.std }, json_file, sort_keys = True, indent = 4)

    def plot_reference_distribution(self, counts):
        visualizer.histogram(counts, name = 'kmer count distribution', path = self.get_current_job_directory(), x_label = 'count', y_label = 'number of kmers')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ChromosomeGcContentEstimationJob(map_reduce.GenomeDependentJob):

    _name = 'ChromosomeGcContentEstimationJob'
    _category = 'programming'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = ChromosomeGcContentEstimationJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.gc = {}
        for i in range(0, 101):
            self.gc[i] = {}
        if c.resume_from_reduce:
            return
        self.chromosomes = extract_whole_genome()
        self.load_reference_counts_provider() 
        #self.exons = bed.load_tracks_from_file_as_dict('/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/gc_profile.bed.Clean')
        self.round_robin(self.chromosomes)

    def transform(self, sequence, chrom):
        l = len(sequence)
        d = l / 100
        t = time.time()
        for i in range(0, d):
            window = sequence[i * 100: (i + 1) * 100]
            gc = calculate_gc_content(window)
            #print(window, gc)
            kmers = c_extract_kmers(32, self.reference_counts_provider.get_kmer_count, 1, True, True, window)
            for kmer in kmers:
                self.gc[gc][kmer] = 0
                break
            if i % 100 == 0:
                s = time.time()
                p = (len(sequence) - i * 100) / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
        return None

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.gc, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def reduce(self):
        c = config.Configuration()
        self.num_threads = 24
        self.kmers = {}
        random.seed(c.seed)
        for batch in self.load_output():
            for i in range(0, 101):
                for kmer in batch[str(i)]:
                    self.gc[i][kmer] = 0
        for gc in self.gc:
            print(gc, len(self.gc[gc]))
            if len(self.gc[gc]) > 10000:
                k = list(self.gc[gc].keys())
                for kmer in [k[i] for i in sorted(random.sample(xrange(len(k)), 10000))]:
                    self.kmers[kmer] = 0
            else:
                for kmer in self.gc[gc]:
                    self.kmers[kmer] = 0
        print('counting', len(self.kmers), 'kmers')
        job = ChromosomeGcContentEstimationJob.CounterHelper(resume_from_reduce = True)
        job.kmers = self.kmers
        kmers = job.execute()
        for gc in self.gc:
            for kmer in self.gc[gc]:
                if kmer in kmers:
                    self.gc[gc][kmer] = kmers[kmer]
        for gc in self.gc:
            for kmer in self.gc[gc]:
                if kmer in kmers:
                    self.gc[gc][kmer] = kmers[kmer]
        # now we have the counts
        self.coverage = []
        for gc in range(0, 96):
            counts = [self.gc[gc][kmer] for kmer in list(filter(lambda k: self.gc[gc][k], self.gc[gc]))]
            mean = numpy.mean(counts)
            std = numpy.std(counts)
            #
            #self.counts = list(filter(lambda x: x < 3 * mean, counts))
            #mean = numpy.mean(counts)
            #std = numpy.std(counts)
            #
            self.counts = list(filter(lambda x: x < 3 * mean, counts))
            mean = numpy.mean(counts)
            std = numpy.std(counts)
            self.coverage.append(mean)
            print(gc, mean)
        for i in range(1, 95):
            if self.coverage[i] > (self.coverage[i - 1] + self.coverage[i + 1]) / 2.0 + 5:
                self.coverage[i] = (self.coverage[i - 1] + self.coverage[i + 1]) / 2.0 + 5
        with open(os.path.join(self.get_current_job_directory(), 'coverage.json'), 'w') as json_file:
            json.dump(self.coverage, json_file)
        visualizer.bar(list(range(0, 96)), [self.coverage], 'gc_corrected_coverage', self.get_current_job_directory(), 'GC Percentage', 'Coverage')

    class CounterHelper(map_reduce.GenomeDependentJob, counter.BaseExactCountingJob):

        _name = 'ChromosomeGcContentEstimationJob'
        _category = 'programming'
        _previous_job = None
        _counter_mode = 3

        @staticmethod
        def launch(**kwargs):
            job = ChromosomeGcContentEstimationJob.CounterHelper(**kwargs)
            job.execute()

        def load_inputs(self):
            self.num_threads = 4
            with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4)
            self.round_robin()

        def reduce(self):
            self.merge_accelerator_counts()
            return self.kmers

        def merge_accelerator_counts(self):
            for i in range(0, self.num_threads):
                print('adding batch', i)
                path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
                self.kmers = {}
                with open (path, 'r') as json_file:
                    line = json_file.readline()
                    while line:
                        kmer = line[:line.find(':')]
                        i = line.find(':')
                        j = line.find(':', i + 1)
                        count = int(line[i + 1: j])
                        total = int(line[j + 1:])
                        canon = canonicalize(kmer)
                        if not canon in self.kmers:
                            self.kmers[canon] = 0
                        self.kmers[canon] += count / 2
                        line = json_file.readline()

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
