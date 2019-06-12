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

from nebula import (
    bed,
    config,
    counter,
    map_reduce,
    statistics,
    visualizer,
)

from nebula.debug import *
from nebula.kmers import *
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

import numpy

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class UniqueKmersDepthOfCoverageEstimationJob(map_reduce.GenomeDependentJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'UniqueKmersDepthOfCoverageEstimationJob'
    _category = 'preprocessing'
    _previous_job = None
    _counter_mode = 3

    @staticmethod
    def launch(**kwargs):
        job = UniqueKmersDepthOfCoverageEstimationJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.load_depth_of_coverage_kmers()
        self.export_accelerator_input()
        self.round_robin()

    def load_depth_of_coverage_kmers(self):
        n = 0
        self.acc = {}
        self.kmers = {}
        for i in range(1, 2):
            self.acc[i] = 0
        self.load_reference_counts_provider() 
        for kmer, count in self.reference_counts_provider.stream_kmers():
            if count in self.acc:
                if self.acc[count] == 1000000:
                    self.acc.pop(count, None)
                    if len(self.acc) == 0:
                        break
                    continue
                self.kmers[canonicalize(kmer)] = 0 
                self.acc[count] += 1
                if self.acc[count] % 1000 == 0:
                    print(self.acc[count], 'kmers with a count of', count, 'so far')
        print('Counting', green(len(self.kmers)), 'unique kmers')
        self.unload_reference_counts_provider()

    def export_accelerator_input(self):
        with open(os.path.join(self.get_current_job_directory(), 'pre_inner_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        del self.kmers

    def merge_count(self, kmer, tokens):
        count = tokens[0] 
        canon = canonicalize(kmer)
        if not canon in self.kmers:
            self.kmers[canon] = 0
        self.kmers[canon] += count / 2

    def reduce(self):
        c = config.Configuration()
        self.kmers = {}
        self.merge_counts()
        self.counts = list(map(lambda kmer: self.kmers[kmer], self.kmers))
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
        # filter outliers
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
# ============================================================================================================================ #
# ============================================================================================================================ #

class ChromosomeGcContentEstimationJob(map_reduce.GenomeDependentJob):

    _name = 'ChromosomeGcContentEstimationJob'
    _category = 'preprocessing'
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
        self.window_size = 96
        for i in range(0, self.window_size + 1):
            self.gc[i] = {}
        self.chromosomes = extract_whole_genome()
        self.load_reference_counts_provider() 
        self.round_robin(self.chromosomes)

    def transform(self, sequence, chrom):
        l = len(sequence)
        d = l / self.window_size
        t = time.time()
        for i in range(0, d):
            window = sequence[i * self.window_size: (i + 1) * self.window_size]
            gc = calculate_gc_content(window)
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
        self.kmers = {}
        random.seed(c.seed)
        for batch in self.load_output():
            for i in range(0, self.window_size + 1):
                for kmer in batch[str(i)]:
                    self.gc[i][kmer] = 0
        for gc in self.gc:
            print(len(self.gc[gc]), 'kmers with GC content', gc)
            if len(self.gc[gc]) > 10000:
                k = list(self.gc[gc].keys())
                for kmer in [k[i] for i in sorted(random.sample(xrange(len(k)), 10000))]:
                    self.kmers[kmer] = gc
            else:
                for kmer in self.gc[gc]:
                    self.kmers[kmer] = gc
        print('counting', len(self.kmers), 'kmers')
        with open(ps.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.gc, json_file, indent = 4)
        return self.gc
