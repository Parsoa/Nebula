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
    map_reduce,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce Job to find reads containing kmers from structural variation events that produce too many novel kmers.
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseExactCountingJob(map_reduce.Job):

    # ============================================================================================================================ #
    # job-specific stuff
    # ============================================================================================================================ #

    def parse_fastq(self):
        name = None
        #tracemalloc.start()
        HEADER_LINE = 0
        SEQUENCE_LINE = 1
        THIRD_LINE = 2
        QUALITY_LINE = 3
        state = HEADER_LINE
        # need to skip invalid lines
        line = self.fastq_file.readline().upper().strip()
        # for the very rare occasion that the first byte in a segment is a line feed
        if len(line) == 0:
            line = self.fastq_file.readline().upper().strip()
        ahead = self.fastq_file.readline().upper().strip()
        n = 0
        m = 0
        t = time.time()
        while ahead:
            #print(state, line)
            if state == HEADER_LINE:
                if line[0] == '@' and ahead[0] != '@':
                    if self.fastq_file.tell() >= (self.index + 1) * self.fastq_file_chunk_size:
                        print(self.index, 'reached segment boundary')
                        break
                    name = line[:-1] # ignore the EOL character
                    state = SEQUENCE_LINE
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
                    print('{:2d}'.format(self.index), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
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

    def round_robin(self):
        for i in range(0, self.num_threads):
            self.batch[i] = {}

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
        for read, name in self.parse_fastq():
            kmers = extract_kmers(c.ksize, read)
            for kmer in kmers:
                if kmer in self.kmers: 
                    self.kmers[kmer]['count'] += 1

    def merge_counts(self):
        c = config.Configuration()
        print('merging kmer counts ...')
        kmers = {}
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            if not os.path.isfile(path):
                print(red('couldn\'t find batch'), i, red('results will be unreliable'))
                continue
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    k = find_kmer(kmer, kmers)
                    if k:
                        kmers[k]['count'] += batch[kmer]['count']
                    else:
                        kmers[kmer] = batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)
        return kmers

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce Job to find reads containing kmers from structural variation events that produce too many novel kmers.
# ============================================================================================================================ #
# ============================================================================================================================ #

class SimulationExactCountingJob(BaseExactCountingJob):

    def parse_fastq(self):
        name = None
        #tracemalloc.start()
        HEADER_LINE = 0
        SEQUENCE_LINE = 1
        THIRD_LINE = 2
        QUALITY_LINE = 3
        state = HEADER_LINE
        # need to skip invalid lines
        line = self.fastq_file.readline().upper().strip()
        # for the very rare occasion that the first byte in a segment is a line feed
        n = 0
        m = 0
        t = time.time()
        while line:
            #print(state, line)
            if state == HEADER_LINE:
                name = line[:-1] # ignore the EOL character
                state = SEQUENCE_LINE
            elif state == SEQUENCE_LINE:
                state = THIRD_LINE
                seq = line[:-1] # ignore the EOL character
                n += 1
                if n == 100000:
                    n = 0
                    m += 1
                    c = self.fastq_file.tell()
                    s = time.time()
                    p = c / float(os.path.getsize(self.fastq_file.name))
                    e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                    print('{:2d}'.format(self.index), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
                yield seq, name
            elif state == THIRD_LINE:
                state = QUALITY_LINE
            elif state == QUALITY_LINE:
                state = HEADER_LINE
            line = self.fastq_file.readline()
        print(self.index, ' end of input')

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def round_robin(self):
        chroms = range(1, 23)
        chroms.append('x')
        chroms.append('y')
        chroms = ['chr' + str(chrom) for chrom in chroms]
        files = {}
        for chrom in chroms:
            for i in range (1, 3):
                for j in range(1, 3):
                    path = os.path.join(self.get_simulation_directory(), chrom + '_strand_' + str(i) + '.' + str(j) + '.fq')
                    files[path] = path
                    print(path)
        map_reduce.Job.round_robin(self, files)

    def run_batch(self, batch):
        c = config.Configuration()
        for track in batch:
            try:
                self.transform(batch[track], track)
            except Exception as e:
                print(red(e))
                traceback.print_exc()
        self.output_batch(self.kmers)

    def transform(self, track, track_name):
        c = config.Configuration()
        self.fastq_file = open(track, 'r')
        for read, name in self.parse_fastq():
            kmers = extract_kmers(read)
            for kmer in kmers:
                if kmer in self.kmers: 
                    self.kmers[kmer]['count'] += 1

