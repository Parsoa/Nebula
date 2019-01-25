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

class UniqueKmersDepthOfCoverageEstimationJob(map_reduce.GenomeDependentJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'UniqueKmersDepthOfCoverageEstimationJob'
    _category = 'programming'
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
        with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        self.round_robin()

    def merge_count(self, kmer, tokens):
        count = tokens[0] 
        canon = canonicalize(kmer)
        self.kmers[canon]['count'] += count / 2

    def reduce(self):
        c = config.Configuration()
        self.merge_counts()
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
# MapReduce job that extracts the kmers for the intervals specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class GappedKmersDepthOfCoverageEstimationJob(map_reduce.GenomeDependentJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'GappedKmersDepthOfCoverageEstimationJob'
    _category = 'programming'
    _previous_job = None
    _counter_mode = 1

    @staticmethod
    def launch(**kwargs):
        job = GappedKmersDepthOfCoverageEstimationJob(**kwargs)
        job.execute()

    def pre_process(self):
        #job = GappedKmersDepthOfCoverageEstimationJob.ExtractExonGappedKmersHelper()
        #job.execute()
        #job = GappedKmersDepthOfCoverageEstimationJob.ScoreExonGappedKmersHelper()
        #job.execute()
        pass

    # ============================================================================================================================ #
    # ============================================================================================================================ #
    # Extract Kmers
    # ============================================================================================================================ #
    # ============================================================================================================================ #

    class ExtractExonGappedKmersHelper(map_reduce.GenomeDependentJob):

        _name = 'GappedKmersDepthOfCoverageEstimationJob'
        _category = 'programming'
        _previous_job = None

        def load_inputs(self):
            c = config.Configuration()
            exons = bed.load_tracks_from_file_as_dict(c.exons)
            self.round_robin(exons)
            extract_whole_genome()
            self.kmers = {}

        def transform(self, track, track_name):
            chromosome = extract_chromosome(track.chrom.lower())
            sequence = chromosome[track.begin: track.end]
            for i in range(0, len(sequence) - 50):
                kmer = canonicalize(sequence[i: i + 37])
                self.kmers[kmer] = 0
            return None

        def output_batch(self, batch):
            json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
            json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
            json_file.close()
            exit()

    # ============================================================================================================================ #
    # ============================================================================================================================ #
    # Scoring Kmers
    # ============================================================================================================================ #
    # ============================================================================================================================ #

    class ScoreExonGappedKmersHelper(map_reduce.GenomeDependentJob):

        _name = 'GappedKmersDepthOfCoverageEstimationJob'
        _category = 'programming'
        _previous_job = None

        def load_inputs(self):
            c = config.Configuration()
            self.kmers = {}
            self.half_mers = {}
            with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers:
                    left = kmer[:c.hsize]
                    right = kmer[-c.hsize:]
                    self.half_mers[left] = True
                    #self.half_mers[right] = True
                    #self.half_mers[reverse_complement(left)] = True
                    self.half_mers[reverse_complement(right)] = True
                    self.kmers[left + right] = 0
            print('scoring', len(self.kmers), 'kmers')
            print('scoring', len(self.half_mers), 'kmers')
            self.chroms = extract_whole_genome()
            self.round_robin(self.chroms)

        def transform(self, sequence, chrom):
            c = config.Configuration()
            t = time.time()
            l = len(sequence)
            for index in range(0, l - 50):
                if index % 100000 == 1:
                    s = time.time()
                    p = index / float(len(sequence))
                    e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                    print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
                half_mer = sequence[index: index + c.hsize]
                if not half_mer in self.half_mers:
                    continue
                for k in range(32, 38):
                    kmer = canonicalize(sequence[index: index + k])
                    kmer = kmer[:c.hsize] + kmer[-c.hsize:]
                    if kmer in self.kmers:
                        self.kmers[kmer] += 1
    
        def output_batch(self, batch):
            json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
            json.dump(self.kmers, json_file, sort_keys = True, indent = 4)
            json_file.close()
            exit()

        def merge_counts(self):
            c = config.Configuration()
            print('merging kmer counts ...')
            for i in range(0, self.num_threads):
                print('adding batch', i)
                path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
                with open (path, 'r') as json_file:
                    batch = json.load(json_file)
                    for kmer in batch:
                        self.kmers[kmer] += batch[kmer]

        def reduce(self):
            self.merge_counts()
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4)

    # ============================================================================================================================ #
    # ============================================================================================================================ #
    # Main job, count kmers in sample
    # ============================================================================================================================ #
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = {}
        self.half_mers = {}
        path = os.path.join(self.get_current_job_directory(), 'gapped_kmers.json')
        print(path)
        with open(path, 'r') as json_file:
            kmers = json.load(json_file)
        kmers = list(filter(lambda k: kmers[k] == 1, kmers))
        print(len(kmers))
        if len(kmers) > 100000:
            kmers = [kmers[i] for i in sorted(random.sample(xrange(len(kmers)), 100000))]
        print('Counting', cyan(len(kmers)), 'gapped kmers')
        for kmer in kmers:
            left = kmer[:c.hsize]
            right = kmer[-c.hsize:]
            self.kmers[kmer] = {"gap": 5} 
            if not left in self.half_mers:
                self.half_mers[left] = {}
            self.half_mers[left][right] = kmer
            left = reverse_complement(left)
            right = reverse_complement(right)
            if not right in self.half_mers:
                self.half_mers[right] = {}
            self.half_mers[right][left] = kmer
        self.export_accelerator_input()
        self.round_robin()

    def export_accelerator_input(self):
        with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
            json.dump(self.half_mers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)

    def merge_count(self, kmer, tokens):
        count = tokens[0]
        if not kmer in self.kmers:
            self.kmers[kmer] = 0
        self.kmers[kmer] += count

    def reduce(self):
        c = config.Configuration()
        self.kmers = {}
        self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4, sort_keys = True)
        #
        counts = [self.kmers[kmer] for kmer in self.kmers] 
        self.mean = numpy.mean(counts)
        self.std = numpy.std(counts)
        print(len(counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter outliers
        self.counts = list(filter(lambda x: x < 3 * self.mean, counts))
        self.mean = numpy.mean(counts)
        self.std = numpy.std(counts)
        print(len(counts))
        print('mean:', self.mean)
        print('std:', self.std)
        # filter outliers
        self.counts = list(filter(lambda x: x < 2 * self.mean, counts))
        self.mean = numpy.mean(counts)
        self.std = numpy.std(counts)
        print(len(counts))
        print('mean:', self.mean)
        print('std:', self.std)

# ============================================================================================================================ #
# ============================================================================================================================ #
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
        self.round_robin(self.chromosomes)

    def transform(self, sequence, chrom):
        l = len(sequence)
        d = l / 100
        t = time.time()
        for i in range(0, d):
            window = sequence[i * 100: (i + 1) * 100]
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
        with open(os.path.join(self.get_current_job_directory(), 'gc.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4)
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
            #
            self.coverage.append(mean)
            print(gc, mean)
        for i in range(1, 95):
            if self.coverage[i] > (self.coverage[i - 1] + self.coverage[i + 1]) / 2.0 + 5:
                self.coverage[i] = (self.coverage[i - 1] + self.coverage[i + 1]) / 2.0 + 5
        with open(os.path.join(self.get_current_job_directory(), 'coverage.json'), 'w') as json_file:
            json.dump(self.coverage, json_file)
        visualizer.bar(list(range(0, 96)), [self.coverage], 'gc_corrected_coverage', self.get_current_job_directory(), 'GC Percentage', 'Coverage')

    # ============================================================================================================================ #
    # ============================================================================================================================ #
    # Scoring Kmers
    # ============================================================================================================================ #
    # ============================================================================================================================ #

    class CounterHelper(map_reduce.GenomeDependentJob, counter.BaseExactCountingJob):

        _name = 'ChromosomeGcContentEstimationJob'
        _category = 'programming'
        _previous_job = None
        _counter_mode = 3

        def load_inputs(self):
            self.num_threads = 4
            with open(os.path.join(self.get_current_job_directory(), 'pre_kmers.json'), 'w') as json_file:
                json.dump(self.kmers, json_file, indent = 4)
            self.round_robin()
            self.num_threads = 4

        def merge_count(self, kmer, tokens):
            canon = canonicalize(kmer)
            count = tokens[0]
            # randomness in selection earlier
            if not canon in self.kmers:
                self.kmers[canon] = 0
            self.kmers[canon] += count / 2

        def reduce(self):
            self.kmers = {}
            self.merge_counts()
            return self.kmers

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
