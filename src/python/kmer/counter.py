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

class RecursiveMergeCountsJob(map_reduce.Job):

    @staticmethod
    def launch(previous_job_directory, **kwargs):
        job = RecursiveMergeCountsJob(job_name = 'RecursiveMergeCountsJob_', previous_job_name = previous_job_directory.split('/')[-1] + '_', previous_job_directory = previous_job_directory, **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # Helper to speed up exporting
    # ============================================================================================================================ #

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        pass

    def distribute_workload(self):
        for index in range(0, self.num_threads, self.batch_size):
            pid = os.fork()
            if pid == 0:
                # forked process
                self.index = index
                self.transform()
                exit()
            else:
                # main process
                self.children[pid] = index
                print('spawned child', '{:2d}'.format(index), ':', pid)
        print(cyan('done distributing workload'))

    def transform(self):
        c = config.Configuration()
        base = self.index
        child = int(self.index + self.batch_size / 2)
        print('merging ', blue(base), 'and', green(child))
        base_dir = self.get_previous_job_directory() if self.batch_size == 1 else self.get_current_job_directory()
        suffix = '' if self.batch_size == 1 else '_' + str(self.batch_size / 2)
        with open(os.path.join(base_path, 'batch_' + str(base) + suffix + '.json'), 'r') as json_file:
            kmers = json.load(json_file)
        with open(os.path.join(base_path, 'batch_' + str(child) + suffix + '.json'), 'r') as json_file:
            batch = json.load(json_file)
            for kmer in batch:
                kmers[kmer] += batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'batch_' + str(base) + '_' + str(self.batch_size) + '.json'), 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4)

    def reduce(self):
        pass

    def get_output_directory(self):
        return self.previous_job_directory

    def get_previous_job_directory(self):
        return self.previous_job_directory

    def get_current_job_directory(self):
        # get rid of the final _
        return os.path.abspath(os.path.join(self.get_output_directory(), self.job_name[:-1]))

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

    #def merge_counts_recursive(self):
    #    for i in range(1, math.ceil(math.log(self.num_threads, 2)) + 1):
    #        job = RecursiveMergeCountsJob.launch(batch_size = 2 ** i, previous_job_directory = self.get_current_job_directory())

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

