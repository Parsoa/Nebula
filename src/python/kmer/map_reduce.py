import io
import os
import re
import pwd
import sys
import copy
import math
import time
import atexit
import argparse
import traceback
import tracemalloc

from kmer import (
    config,
)

from kmer.kmers import *
from kmer.commons import *
from kmer.commons import pretty_print as print

import colorama
import memory_profiler

import rapidjson as json

def on_exit(job):
    print(colorama.Fore.GREEN + 'job', job.index, 'exiting', colorama.Fore.WHITE)

print('importing map_reduce.py')
# ============================================================================================================================ #
# Job class, describes a MapReducer job
# Will apply a transformation to a library of structural variation:
# 1. Divides the library into a number of `batches`
# 2. Each batch includes one or more `tracks` (each track is a structural variation)
# 3. Applies the function `transform` to each track and outputs the result
# 4. Merges the transformed tracks into a whole
# ============================================================================================================================ #

class Job(object):

    def __init__(self, job_name, previous_job_name, **kwargs):
        self.job_name = job_name
        self.previous_job_name = previous_job_name
        self.index = -1
        self.batch = {}
        self.children = {}
        self.run_for_certain_batches_only = False
        self.resume_from_reduce = False
        for k, v in kwargs.items():
            setattr(self, k, v)

    def prepare(self):
        pass

    def check_cli_arguments(self, args):
        pass

    def execute(self):
        c = config.Configuration()
        self.check_cli_arguments(None)
        self.prepare()
        self.create_output_directories()
        self.find_thread_count()
        if not self.resume_from_reduce:
            print('normal execution flow')
            self.load_inputs()
            self.distribute_workload()
            self.wait_for_children()
        else:
            print('resuming from reduce')
        self.reduce()
        # self.clean_up()

    # this for when you need to make small adjustments to the output after the job has finished but don't want to run it all over again
    def post_process(self):
        pass

    # this is upperbounded by --threads cli argument
    def find_thread_count(self):
        c = config.Configuration()
        max_index = 0
        for index in range(0, c.max_threads):
            path = os.path.join(self.get_previous_job_directory(), 'batch_' + str(index) + '.json')
            if os.path.isfile(path):
                max_index = index + 1
        self.num_threads = max_index

    def load_inputs(self):
        for index in range(0, self.num_threads):
            path = os.path.join(self.get_previous_job_directory(), 'batch_' + str(index) + '.json')
            with open(path, 'r') as json_file:
                self.batch[index] = json.load(json_file)

    def distribute_workload(self):
        for index in range(0, self.num_threads):
            if self.run_for_certain_batches_only:
                if not index in self.batches_to_run:
                    continue
            pid = os.fork()
            if pid == 0:
                # forked process
                self.index = index
                # atexit.register(on_exit, self)
                self.run_batch(self.batch[index])
                exit()
            else:
                # main process
                self.children[pid] = index
                print('spawned child', '{:2d}'.format(index), ':', pid)
        print(colorama.Fore.CYAN + 'done distributing workload')

    def run_batch(self, batch):
        c = config.Configuration()
        remove = {}
        n = 0
        start = time.time()
        for track in batch:
            batch[track] = self.transform(batch[track], track)
            if batch[track] == None:
                remove[track] = True
            n = n + 1
            t = time.time()
            p = float(n) / len(batch)
            eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
            print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
        for track in remove:
            batch.pop(track, None)
        # ths forked process will exit after the following function call
        self.output_batch(batch)

    def transform(self, track, track_name):
        return track

    # This MUST call exit()
    def output_batch(self, batch):
        # output manually, io redirection could get entangled with multiple client/servers
        n = 0
        while False:
            if self.index == 0:
                break
            if os.path.isfile(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index - 1) + '.json')):
                print('found output for', self.index - 1)
                break
            n += 1
            if n == 100000:
                print(self.index, 'waiting for', self.index - 1)
                n = 0 
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(batch, json_file, sort_keys = True, indent = 4)
        json_file.close()
        exit()

    def wait_for_children(self):
        while True:
            if len(self.children) == 0:
                break
            (pid, e) = os.wait()
            index = self.children[pid]
            self.children.pop(pid, None)
            if os.path.isfile(os.path.join(self.get_current_job_directory(), 'batch_' + str(index) + '.json')):
                print(red('pid', '{:5d}'.format(pid), ', index', '{:2d}'.format(index), 'finished,', '{:2d}'.format(len(self.children)), 'remaining'))
            else:
                print(red('pid', '{:5d}'.format(pid), ', index', '{:2d}'.format(index), 'finished didn\'t produce output,', len(self.children), 'remaining'))
        print(cyan('all forks done, merging output ...'))

    def plot(self, outputs):
        pass

    def sort(self, outputs):
        pass

    def reduce(self):
        c = config.Configuration()
        output = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    output.update(batch)
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)
        self.sort(output)
        self.plot(output)

    def clean_up(self):
        c = config.Configuration()
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            os.remove(path)

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/' + str(c.ksize) + '/'))

    def get_previous_job_directory(self):
        # get rid of the final _
        return os.path.abspath(os.path.join(self.get_output_directory(), self.previous_job_name[:-1]))

    def get_current_job_directory(self):
        # get rid of the final _
        return os.path.abspath(os.path.join(self.get_output_directory(), self.job_name[:-1]))

    def create_output_directories(self):
        dir = self.get_output_directory()
        if not os.path.exists(dir):
            os.makedirs(dir)
        dir = self.get_current_job_directory()
        if not os.path.exists(dir):
            os.makedirs(dir)

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce Job to find reads containing kmers from structural variation events that produce too many novel kmers.
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseExactCountingJob(Job):

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
        line = self.fastq_file.readline().strip()
        ahead = self.fastq_file.readline().strip()
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
                    #snapshot = tracemalloc.take_snapshot()
                    #display_top(snapshot)
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

    def check_cli_arguments(self, args):
        #tracemalloc.start()
        # --bed to specify the set of structural variations
        # --fastq: the genome from which we are getting the kmer counts
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

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
                    self.kmers[kmer] += 1

    def merge_counts(self):
        c = config.Configuration()
        # merge kmer counts from all children
        print('merging kmer counts ...')
        kmers = {}
        for i in range(0, self.num_threads):
            print('batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json') 
            if not os.path.isfile(path):
                print(colorama.Fore.RED + 'couldn\'t find batch', i, ' results will be suspicious')
                continue
            with open (path, 'r') as json_file:
                batch = json.load(json_file)
                for kmer in batch:
                    if not kmer in kmers:
                        kmers[kmer] = 0
                    kmers[kmer] += batch[kmer]
        return kmers
