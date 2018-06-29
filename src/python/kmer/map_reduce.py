from __future__ import print_function

import io
import gc
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import atexit
import psutil
import argparse
import traceback

from shutil import copyfile

from kmer import (
    config,
    counttable,
)

from kmer.sv import StructuralVariation, Inversion, Deletion, SNP
from kmer.kmers import *
from kmer.commons import *
print = pretty_print

def on_exit(job):
    print(green('job', job.index, 'exiting'))

# ============================================================================================================================ #
# Job class, describes a MapReducer job
# Will apply a transformation to a library of structural variation:
# 1. Divides the library into a number of `batches`
# 2. Each batch includes one or more `tracks` (each track is a structural variation)
# 3. Applies the function `transform` to each track and outputs the result
# 4. Merges the transformed tracks into a whole
# ============================================================================================================================ #

class Job(object):

    def __init__(self, job_name, previous_job_name, category = 'output', **kwargs):
        self.job_name = job_name
        self.previous_job_name = previous_job_name
        self.category = category
        self.index = -1
        self.batch = {}
        self.batch_file_prefix = 'batch_'
        self.previous_job_batch_file_prefix = 'batch_'
        self.children = {}
        self.run_for_certain_batches_only = False
        self.resume_from_reduce = False
        for k, v in kwargs.items():
            setattr(self, k, v)

    def prepare(self):
        # Declare and initialize needed variables here
        pass

    def check_cli_arguments(self, args):
        # Check if every needed argument is passed and in good form
        pass

    def execute(self):
        c = config.Configuration()
        self.check_cli_arguments(None)
        self.create_output_directories()
        self.find_thread_count()
        if not self.resume_from_reduce:
            print('normal execution flow')
            self.prepare()
            self.load_inputs()
            self.distribute_workload()
            self.wait_for_children()
        else:
            print('resuming from reduce')
        output = self.reduce()
        self.plot(output)

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def prepare(self):
        pass

    def load_inputs(self):
        tracks = self.load_previous_job_results()
        self.round_robin(tracks)

    def load_previous_job_results(self):
        path = os.path.join(self.get_previous_job_directory(), self.previous_job_batch_file_prefix + 'merge.json')
        print(path)
        with open(path, 'r') as json_file:
            return json.load(json_file)

    def round_robin(self, tracks, name_func = lambda x: x, filter_func = lambda x: False):
        c = config.Configuration()
        n = 0
        for track in tracks:
            track_name = name_func(track)
            if filter_func(tracks[track]):
                continue
            index = n % c.max_threads 
            if not index in self.batch:
                self.batch[index] = {}
            self.batch[index][track_name] = tracks[track]
            print(blue('assigned ', track_name, ' to ', index))
            n = n + 1
            self.num_threads = min(c.max_threads, n)

    def distribute_workload(self):
        for index in range(0, self.num_threads):
            if self.run_for_certain_batches_only:
                if not index in self.batches_to_run:
                    continue
            pid = os.fork()
            if pid == 0:
                self.index = index
                # atexit.register(on_exit, self)
                self.run_batch(self.batch[index])
                exit()
            else:
                self.children[pid] = index
                print('spawned child', '{:2d}'.format(index), ':', pid)
        print(cyan('done distributing workload'))

    def run_batch(self, batch):
        c = config.Configuration()
        remove = {}
        n = 0
        start = time.time()
        for track in batch:
            try:
                batch[track] = self.transform(batch[track], track)
                if batch[track] == None:
                    remove[track] = True
            except Exception as e:
                print(red(e))
                traceback.print_exc()
                remove[track] = True
            n = n + 1
            t = time.time()
            p = float(n) / len(batch)
            eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
            print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
            gc.collect()
        for track in remove:
            batch.pop(track, None)
        # if there is no output, don't write anything
        if not batch:
            exit()
        self.output_batch(batch)

    def transform(self, track, track_name):
        return track

    # This MUST call exit()
    def output_batch(self, batch):
        n = 0
        json_file = open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + str(self.index) + '.json'), 'w')
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
            if os.path.isfile(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + str(index) + '.json')):
                print(red('pid', '{:5d}'.format(pid) + ', index', '{:2d}'.format(index), 'finished,', '{:2d}'.format(len(self.children)), 'remaining'))
            else:
                print(red('pid', '{:5d}'.format(pid) + ', index', '{:2d}'.format(index), 'finished didn\'t produce output,', len(self.children), 'remaining'))
        print(cyan('all forks done, merging output ...'))

    def reduce(self):
        c = config.Configuration()
        output = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    output.update(batch)
        with open(os.path.join(self.get_current_job_directory(), self.batch_file_prefix + 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)
        return output

    def load_output_batch(self, index):
        path = os.path.join(self.get_current_job_directory(), self.batch_file_prefix + str(index) + '.json')
        with open(path, 'r') as json_file:
            output = json.load(json_file)
            return output

    def plot(self, output):
        pass

    def sort(self, output):
        pass

    def post_process(self):
        pass

    def clean_up(self):
        c = config.Configuration()
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            os.remove(path)

    # ============================================================================================================================ #
    # misc helpers
    # ============================================================================================================================ #

    def get_sv_type(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        if bed_file_name.find('DEL') != -1:
            return Deletion
        if bed_file_name.find('INV') != -1:
            return Inversion
        return Deletion

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        d = 'simulation' if c.simulation else self.category
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../' + d + '/' + bed_file_name + '/' + str(c.ksize) + '/'))

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
# Base class for every job that is a direct part of the genotyping process
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseGenotypingJob(Job ):

    def get_output_directory(self):
        c = config.Configuration()
        if c.simulation:
            return Job.get_output_directory(self)
        else:
            fastq_file_name = c.fastq_file.split('/')[-1][::-1].split('.')[-1][::-1]
            return os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../' + self.category + '/genotyping/' + fastq_file_name))

    def get_current_job_directory(self):
        c = config.Configuration()
        if c.simulation:
            return Job.get_current_job_directory(self)
        else:
            bed_file_name = c.bed_file.split('/')[-1]
            return os.path.abspath(os.path.join(self.get_output_directory(), self.job_name[:-1], bed_file_name))

    def get_previous_job_directory(self):
        c = config.Configuration()
        print(c.simulation)
        if c.simulation:
            return Job.get_previous_job_directory(self)
        else:
            bed_file_name = c.bed_file.split('/')[-1]
            return os.path.abspath(os.path.join(self.get_output_directory(), self.previous_job_name[:-1], bed_file_name))

