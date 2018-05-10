from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import argparse
import operator
import traceback

from kmer import (
    bed,
    sets,
    config,
    counttable,
    map_reduce,
    statistics,
)

from kmer.sv import StructuralVariation, Inversion, Deletion, SNP
from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import pybedtools

import plotly.offline as plotly
import plotly.graph_objs as graph_objs

# ============================================================================================================================ #
# ============================================================================================================================ #
# base class for all genotyping jobs
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseTrainingJob(map_reduce.Job):

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../training/' + bed_file_name + '/'))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class TrainingJob(BaseTrainingJob):

    def find_thread_count(self):
        c = config.Configuration()
        return c.max_threads

    def load_input(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'merge.json'), 'r') as json_file:
            tracks = json.load(json_file)
        n = 0
        for track in tracks:
            index = n % c.max_threads 
            if not index in self.batch:
                self.batch[index] = {}
            self.batch[index][track] = tracks[track]
            print(blue('assigned', track, 'to', index))
            n = n + 1
            self.num_threads = max(self.num_threads, index + 1)

    def transform(self, track, track_name):
        seq = simulator.extract_chromosome(track.chome)
        seq = seq[max(0, track.start - 250) : min(len(seq), track.end + 250)]
        kmers = get_canonical_kmers(seq)

    def get_current_job_directory(self):
        c = config.Configuration()
        return os.path.abspath(os.path.join(self.get_output_directory(), 'TrainingJob'))

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/31/MostLikelyBreakPointsJob/'))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
