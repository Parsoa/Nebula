import io
import os
import re
import pwd
import sys
import copy
import json
import time
import atexit
import argparse
import traceback

from kmer import (
    config,
    commons,
)

import colorama
import memory_profiler

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
        if 'batches_to_run' in kwargs:
            self.run_for_certain_batches_only = True
            self.batches_to_run = kwargs['batches_to_run']
        if 'resume_from_reduce' in kwargs:
            print('resuming from reduce')
            self.resume_from_reduce = True

    def prepare(self):
        pass

    def check_cli_arguments(self, args):
        pass

    def execute(self):
        c = config.Configuration()
        self.prepare()
        self.create_output_directories()
        self.find_thread_count()
        if not self.resume_from_reduce:
            self.load_inputs()
            self.distribute_workload()
            self.wait_for_children()
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
                atexit.register(on_exit, self)
                self.run_batch(self.batch[index])
            else:
                # main process
                self.children[pid] = index
                print('spawned child ', index, ':', pid)

    def run_batch(self, batch):
        c = config.Configuration()
        remove = {}
        for track in batch:
            batch[track] = self.transform(batch[track], track)
            if batch[track] == None:
                remove[track] = True
        for track in remove:
            batch.pop(track, None)
        # ths forked process will exit after the following function call
        self.output_batch(batch)

    def transform(self, track, track_name):
        return track

    @memory_profiler.profile()
    def output_batch_profile(self, batch):
        print('outputting batch with profiling', self.index)
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
        print('output', self.index, ':', len(batch))
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        #json.dump(batch, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        for chunk in json.JSONEncoder().iterencode(batch):
            json_file.write(chunk)
        json_file.close()
        exit()

    def output_batch(self, batch):
        print('outputting batch', self.index)
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
        print('output', self.index, ':', len(batch))
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(batch, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        #for chunk in json.JSONEncoder().iterencode(batch):
        #    json_file.write(chunk)
        json_file.close()
        exit()

    def wait_for_children(self):
        while True:
            (pid, e) = os.wait()
            index = self.children[pid]
            self.children.pop(pid, None)
            if os.path.isfile(os.path.join(self.get_current_job_directory(), 'batch_' + str(index) + '.json')):
                print(colorama.Fore.RED + 'pid: ', pid, index, 'finished,', len(self.children), 'remaining', colorama.Fore.WHITE)
            else:
                print(colorama.Fore.RED + 'pid: ', pid, index, 'finished didn\'t produce output,', len(self.children), 'remaining', colorama.Fore.WHITE)
            if len(self.children) == 0:
                break
        print('all forks done, merging output ...')

    #TODO: depracate this
    def merge(self, outputs):
        pass

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
            json.dump(output, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        self.merge(output)
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
