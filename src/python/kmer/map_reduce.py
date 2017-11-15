import io
import os
import re
import pwd
import sys
import copy
import json
import time
import argparse
import traceback

# ============================================================================================================================ #
# Job class, describes a MapReducer job
# ============================================================================================================================ #

class Job(object):

    def __init__(self, job_name, previous_job_name, **kwargs):
        self.job_name = job_name
        self.previous_job_name = previous_job_name
        self.index = -1
        self.batch = {}
        self.children = {}

    def prepare():
        pass

    def execute(self):
        c = config.Configuration()
        self.prepare()
        self.find_thread_count()
        self.load_inputs()
        self.distribute_workload()
        self.wait_for_children(children)
        self._reduce()
        self.clean_up()

    def find_thread_count(self):
        c = config.Configuration()
        max_index = 0
        for index in range(0, c.max_threads):
            path = os.path.join(self.get_output_directory(), 'batch_' + self.previous_job_name + str(index) + '.json')
            if os.path.isfile(path):
                max_index = index + 1
        self.num_threads = max_index

    def load_inputs():
        for index in range(0, self.num_threads):
            path = os.path.join(self.get_output_directory(), 'batch_' + self.previous_job_name + str(index) + '.json')
            with open(path, 'r') as json_file:
                self.batch[index] = json.load(json_file)

    def distribute_workload(self):
        for index in range(0, self.num_threads):
            pid = os.fork()
            if pid == 0:
                # forked process
                self.run_batch(self.batch[i])
            else:
                # main process
                self.children[pid] = True
                print('spawned child ', pid)

    def run_batch(self, batch):
        c = config.Configuration()
        for track in batch:
            batch[track] = self.transform(batch[track], track)
        self.output_batch(batch)
        print(colorama.Fore.GREEN, 'process ', self.index, ' done')

    def transform(track, track_name, index):
        return track

    def output_batch(self, batch):
        # output manually, io redirection could get entangled with multiple client/servers
        with open(os.path.join(self.get_output_directory(), 'batch_' + self.job_name + str(index) + '.json'), 'w') as json_file:
            json.dump(batch, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        exit()

    def wait_for_children(self):
        while True:
            (pid, e) = os.wait()
            self.children.pop(pid, None)
            print(colorama.Fore.RED, 'pid ', pid, 'finished')
            if len(self.children) == 0:
                break
        print('all forks done, merging output ...')

    def merge(self, outputs):
        pass

    def plot(self, outputs):
        pass

    def _reduce(self):
        c = config.Configuration()
        output = {}
        for i in range(0, self.num_threads):
            with open(os.path.join(self.get_output_directory(), 'batch_' + self.job_name + str(i) + '.json'), 'r') as json_file:
                batch = json.load(json_file)
                output.update(batch)
        self.merge(output)
        self.plot(ouput)
        with open(os.path.join(self.get_output_directory(), 'merge_' + self.job_name + '.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4, separators = (',', ': '))

    def clean_up(self):
        c = config.Configuration()
        for i in range(0, self.num_threads):
            # might fail because there weren't as many as i processes
            path = os.path.join(self.get_output_directory(), 'batch_' + self.job_name + str(i) + '.json')
            os.remove(path)

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/' + str(c.ksize)))