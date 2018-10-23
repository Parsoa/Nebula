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
    bed,
    config,
    counttable,
)

from kmer.sv import StructuralVariation, Inversion, Deletion
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

    def __init__(self, **kwargs):
        self.index = -1
        self.batch = {}
        self.children = {}
        self.run_for_certain_batches_only = False
        self.resume_from_reduce = False
        for k, v in kwargs.items():
            print('adding attr', green(k), blue(v))
            setattr(self, k, v)

    def prepare(self):
        # Declare and initialize needed variables here before the mpa reduce flow beings
        pass

    def check_cli_arguments(self, args):
        # Check if every needed argument is passed and in good form
        pass

    def execute(self):
        c = config.Configuration()
        self.check_cli_arguments(None)
        self.prepare()
        self.create_output_directories()
        self.find_thread_count()
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

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def prepare(self):
        pass

    def load_inputs(self):
        tracks = self.load_previous_job_results()
        self.round_robin(tracks)

    def load_previous_job_results(self):
        path = os.path.join(self.get_previous_job_directory(),  'batch_merge.json')
        with open(path, 'r') as json_file:
            return json.load(json_file)

    def round_robin(self, tracks, name_func = lambda x: x, filter_func = lambda x: False):
        c = config.Configuration()
        print('round robin', c.max_threads)
        n = 0
        self.batch = {}
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
        self.output_batch(batch)

    def transform(self, track, track_name):
        return track

    # This MUST call exit()
    def output_batch(self, batch):
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

    def load_tracks(self, name = 'all.bed'):
        c = config.Configuration()
        if c.simulation:
            return {str(track): self.get_sv_type()(track) for track in bed.load_tracks_from_file(os.path.join(self.get_simulation_directory(), name))}
        else:
            return {str(track): self.get_sv_type()(track) for track in bed.load_tracks_from_file(c.bed_file)}

    def get_sv_type(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        if bed_file_name.find('DEL') != -1:
            return Deletion
        if bed_file_name.find('INV') != -1:
            return Inversion
        return Deletion

    def load_reference_counts_provider(self):
        c = config.Configuration()
        if c.simulation:
            path = os.path.join(self.get_simulation_directory(), 'reference_' + str(c.ksize) + '.jf')
        else:
            path = os.path.join(c.jellyfish_base, c.reference, 'mer_counts_' + str(c.ksize) + '.jf')
        print('Reference kmer index:', green(path))
        self.reference_counts_provider = counttable.JellyfishCountsProvider(path)

    def get_gapped_reference_counts_provider(self):
        c = config.Configuration()
        g = c.gap
        k = (c.ksize / 2) * 2
        print(g, k)
        if c.simulation:
            return os.path.join(self.get_simulation_directory(), 'reference_' + str(g + k) + '.jf')
        else:
            return os.path.join(c.jellyfish_base, c.reference, 'mer_counts_' + str(g + k) + '.jf')

    def merge_counts(self, *keywords):
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
                        for keyword in keywords:
                            if keyword in kmers[k]:
                                kmers[k][keyword] += batch[kmer][keyword]
                    else:
                        kmers[kmer] = batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)
        return kmers

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        if c.simulation:
            return os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../simulation/' + bed_file_name + '/' + str(c.description) + '/' + str(c.simulation) + 'x/' + str(c.ksize) + 'k/'))
        else:
            return os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../' + self._category + '/' + bed_file_name + '/' + str(c.ksize) + 'k/'))

    def get_previous_job_directory(self):
        if self._previous_job:
            return os.path.abspath(os.path.join(self.get_output_directory(), self._previous_job._name))
        else:
            return self.get_current_job_directory()

    def get_current_job_directory(self):
        return os.path.abspath(os.path.join(self.get_output_directory(), self._name))

    def get_simulation_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../simulation/' + bed_file_name + '/' + str(c.description) + '/' + str(c.simulation) + 'x/' + str(c.ksize) + 'k/Simulation/'))

    def create_output_directories(self):
        dir = self.get_output_directory()
        print(yellow(dir))
        if not os.path.exists(dir):
            os.makedirs(dir)
        dir = self.get_current_job_directory()
        print('creating directory:', green(dir))
        if not os.path.exists(dir):
            os.makedirs(dir)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Base class for every job that is a direct part of the genotyping process
# ============================================================================================================================ #
# ============================================================================================================================ #

class FirstGenotypingJob(Job):

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        if c.simulation:
            return Job.get_output_directory(self)
        else:
            return os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../' + self._category + '/genotyping/' + c.genome))

    def get_current_job_directory(self):
        c = config.Configuration()
        if c.simulation:
            return Job.get_current_job_directory(self)
        else:
            bed_file_name = c.bed_file.split('/')[-1]
            return os.path.abspath(os.path.join(self.get_output_directory(), bed_file_name, self._name))

    def get_previous_job_directory(self):
        c = config.Configuration()
        d = Job.get_output_directory(self)
        print(d)
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(d, self._previous_job._name))

# ============================================================================================================================ #
# ============================================================================================================================ #
# Base class for every job that is a direct part of the genotyping process
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseGenotypingJob(FirstGenotypingJob):

    def get_previous_job_directory(self):
        c = config.Configuration()
        if c.simulation:
            return Job.get_previous_job_directory(self)
        else:
            bed_file_name = c.bed_file.split('/')[-1]
            return os.path.abspath(os.path.join(self.get_output_directory(), bed_file_name, self._previous_job._name))

