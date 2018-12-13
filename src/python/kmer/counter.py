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
        u = 0
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
                u += 1
                if n == 100000:
                    n = 0
                    m += 1
                    c = self.fastq_file.tell() - self.index * self.fastq_file_chunk_size
                    s = time.time()
                    p = c / float(self.fastq_file_chunk_size)
                    e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                    print('{:2d}'.format(self.index), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
                yield seq, name, u
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
        self.transform()
        self.output_batch(self.kmers)

    def transform(self):
        c = config.Configuration()
        cpp_dir = os.path.join(os.path.dirname(__file__), '../../cpp')
        if c.accelerate:
            command = " " + str(self.index) + " " + self.get_current_job_directory() +  " " + c.fastq_file + " " + str(self.num_threads) + " " + str(self._counter_mode) + " " + ("1" if c.debug else "0")
            print(command)
            if c.debug:
                output = subprocess.call(os.path.join(cpp_dir, "counter_s.out") + command, shell = True)
            else:
                output = subprocess.call(os.path.join(cpp_dir, "counter.out") + command, shell = True)
            exit()
        else:
            for read, name, index in self.parse_fastq():
                self.process_read(read, name, index)

    def process_read(self, read, name):
        c = config.Configuration()
        kmers = extract_kmers(c.ksize, read)
        for kmer in kmers:
            if kmer in self.kmers: 
                self.kmers[kmer]['count'] += kmers[kmer]

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
                        j = line.find(':', i + 1)
                        count = int(line[i + 1: j])
                        total = int(line[j + 1:])
                        canon = canonicalize(kmer)
                        self.kmers[canon]['count'] += count / 2
                        self.kmers[canon]['total'] += total / 2
                        line = json_file.readline()
                        n += 1
                    except Exception as e:
                        print(n, i, line)
                        print(e)
                        traceback.print_exc()
                        debug_breakpoint()
                        line = json_file.readline()
                        n += 1

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce Job to find reads containing kmers from structural variation events that produce too many novel kmers.
# ============================================================================================================================ #
# ============================================================================================================================ #

class SimulationExactCountingJob(BaseExactCountingJob):

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
                for j in range(1, 2):
                    name = chrom + '_strand_' + str(i) + '.' + str(j) + '.fq'
                    path = os.path.join(self.get_simulation_directory(), name)
                    files[name] = path
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

    def transform(self, file_path, file_name):
        c = config.Configuration()
        self.fastq_file = open(file_path, 'r')
        self.fastq_file_chunk_size = os.path.getsize(self.fastq_file.name)
        for read, name in self.parse_fastq():
            self.process_read(read, name)
    
