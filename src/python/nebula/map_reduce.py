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
import argparse
import traceback

from shutil import copyfile

from nebula import (
    bed,
    config,
    counttable,
)

from nebula.kmers import *
from nebula.commons import *
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

    def __init__(self, **kwargs):
        c = config.Configuration()
        self.index = -1
        self.batch = {}
        self.children = {}
        self.resume_from_reduce = c.reduce
        for k, v in kwargs.items():
            setattr(self, k, v)

    def execute(self):
        c = config.Configuration()
        self.pre_process()
        self.find_thread_count()
        self.create_output_directories()
        self.prepare()
        self.load_inputs()
        if not self.resume_from_reduce: 
            print('normal execution flow')
            self.distribute_workload()
            self.wait_for_children()
        else:
            print('resuming from reduce')
        output = self.reduce()
        self.plot(output)
        return output

    def pre_process(self):
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.threads

    def prepare(self):
        pass

    def load_inputs(self):
        tracks = self.load_previous_job_results()
        self.round_robin(tracks)

    def load_previous_job_results(self, job = None):
        if not job:
            path = os.path.join(self.get_previous_job_directory(), 'batch_merge.json')
        else:
            path = os.path.join(os.path.join(self.get_current_job_directory(), '..', job._name), 'batch_merge.json')
        with open(path, 'r') as json_file:
            return json.load(json_file)

    def round_robin(self, tracks, filter_func = lambda x: False):
        c = config.Configuration()
        print('round robin', c.threads)
        n = 0
        self.batch = {}
        for track in tracks:
            if filter_func(tracks[track]):
                continue
            index = n % c.threads
            if not index in self.batch:
                self.batch[index] = {}
            self.batch[index][track] = tracks[track]
            print(blue('assigned ', track, ' to ', index))
            n = n + 1
            self.num_threads = min(c.threads, n)

    def distribute_workload(self):
        for index in range(0, self.num_threads):
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
            if n % 1000 == 0:
                gc.collect()
        for track in remove:
            batch.pop(track, None)
        # if there is no output, don't write anything
        self.output_batch(batch)
        self.on_exit_worker()
        exit()

    def on_exit_worker(self):
        pass

    def transform(self, track, track_name):
        return track

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(batch, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def wait_for_children(self):
        while True:
            if len(self.children) == 0:
                break
            (pid, e) = os.wait()
            index = self.children[pid]
            self.children.pop(pid, None)
            if os.path.isfile(os.path.join(self.get_current_job_directory(), 'batch_' + str(index) + '.json')):
                print(red('pid', '{:5d}'.format(pid) + ', index', '{:2d}'.format(index), 'finished,', '{:2d}'.format(len(self.children)), 'remaining'))
            else:
                print(red('pid', '{:5d}'.format(pid) + ', index', '{:2d}'.format(index), 'finished didn\'t produce output,', len(self.children), 'remaining'))
        print(cyan('all forks done, merging output ...'))

    def reduce(self):
        c = config.Configuration()
        output = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    output.update(batch)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)
        return output

    def load_output(self):
        for i in range(0, self.num_threads):
            print('reading batch', i)
            yield self.load_output_batch(i)

    def load_output_batch(self, index):
        path = os.path.join(self.get_current_job_directory(), 'batch_' + str(index) + '.json')
        if not os.path.isfile(path):
            print(yellow('didn\'t find batch', index))
            return {}
        with open(path, 'r') as json_file:
            output = json.load(json_file)
            return output

    def plot(self, output):
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

    def load_tracks(self, name = 'all.bed'):
        c = config.Configuration()
        if c.simulation:
            return bed.load_tracks_from_file_as_dict(os.path.join(self.get_simulation_directory(), name))
        else:
            return bed.load_tracks_from_file_as_dict(c.bed)

    def get_sv_type(self):
        c = config.Configuration()
        bed_file_name = c.bed.split('/')[-1]
        if bed_file_name.find('DEL') != -1:
            return 'DEL'
        if bed_file_name.find('INV') != -1:
            return 'INV'
        if bed_file_name.find('ALU') != -1:
            return 'ALU'
        return 'DEL'

    def load_reference_counts_provider(self):
        c = config.Configuration()
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish)

    def unload_reference_counts_provider(self):
        del self.reference_counts_provider

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_output_directory(self):
        c = config.Configuration()
        return os.path.abspath(c.workdir)

    def get_current_job_directory(self):
        return os.path.abspath(os.path.join(self.get_output_directory(), self._name))

    def get_previous_job_directory(self):
        c = config.Configuration()
        if self._previous_job:
            return os.path.abspath(os.path.join(self.get_output_directory(), self._previous_job._name))
        else:
            return None

    def create_output_directories(self):
        path = self.get_output_directory()
        if not os.path.exists(path):
            os.makedirs(path)
        path = self.get_current_job_directory()
        if not os.path.exists(path):
            os.makedirs(path)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Base class for every job that is a direct part of the genotyping process
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseGenotypingJob(Job):

    pass

# ============================================================================================================================ #
# ============================================================================================================================ #
# First job in the genotyping pipeline, it may have a predecessor non-genotyping job, so the output directory for the previous
# job may be different
# ============================================================================================================================ #
# ============================================================================================================================ #

class FirstGenotypingJob(BaseGenotypingJob):

    pass

# ============================================================================================================================ #
# ============================================================================================================================ #
# Base class for every job that is a direct part of the genotyping process
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenomeDependentJob(BaseGenotypingJob):

    def get_current_job_directory(self):
        c = config.Configuration()
        if c.simulation:
            s = Job.get_current_job_directory(self)
            print(yellow(s))
            return s
        else:
            return os.path.abspath(os.path.join(self.get_output_directory(), self._name))
