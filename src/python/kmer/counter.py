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
        command = " " + str(self.index) + " " + self.get_current_job_directory() +  " " + c.fastq_file + " " + str(self.num_threads) + " " + str(self._counter_mode) + " " + ("1" if c.debug else "0")
        if c.debug:
            output = subprocess.call(os.path.join(cpp_dir, "counter_s.out") + command, shell = True)
        else:
            output = subprocess.call(os.path.join(cpp_dir, "counter.out") + command, shell = True)
        exit()

    def merge_count(self, kmer, tokens):
        pass

    def merge_counts(self):
        for i in range(0, self.num_threads):
            print('adding batch', i)
            path = os.path.join(self.get_current_job_directory(), 'c_batch_' + str(i) + '.json') 
            with open (path, 'r') as json_file:
                line = json_file.readline()
                while line:
                    tokens = line.split(':')
                    self.merge_count(tokens[0], [int(t) for t in tokens[1:]])
                    line = json_file.readline()
