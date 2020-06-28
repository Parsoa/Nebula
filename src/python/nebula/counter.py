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
    map_reduce,
)

from nebula.kmers import *
from nebula.logger import *

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

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = 1

    def round_robin(self):
        for i in range(0, self.num_threads):
            self.batch[i] = {}

    def run_batch(self, batch):
        c = config.Configuration()
        self.transform()

    def transform(self):
        c = config.Configuration()
        self.create_output_directories()
        cpp_dir = os.path.join(os.path.dirname(__file__), '../../cpp/counter')
        if c.bam:
            command = os.path.join(cpp_dir, "counter.out") + " " + self.get_current_job_directory() +  " " + c.bam + " " + str(16)
        else:
            command = os.path.join(cpp_dir, "counter.out") + " " + self.get_current_job_directory() +  " " + c.fastq[0] + " " + str(16)
        print(command)
        output = subprocess.call(command, shell = True)

    def merge_count(self, kmer, tokens):
        pass

    def merge_counts(self):
        c = config.Configuration()
        with open (os.path.join(self.get_current_job_directory(), 'counts.json'), 'r') as json_file:
            line = json_file.readline()
            while line:
                tokens = line.split(':')
                self.merge_count(tokens[0], [int(t) for t in tokens[1:]])
                line = json_file.readline()
        print('Done aggregating kmer counts.')

