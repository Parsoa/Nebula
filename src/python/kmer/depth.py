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
        #self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
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
        with open(os.path.join(self.get_current_job_directory(), 'stats_' + str(c.ksize) + '.json'), 'w') as json_file:
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
