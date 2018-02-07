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

from multiprocessing import Process

from kmer import (
    bed,
    sets,
    config,
    commons,
    counttable,
    map_reduce,
    statistics,
    count_server,
)

from kmer.kmers import *

from kmer.prune import CountKmersExactJob


# ============================================================================================================================ #
# ============================================================================================================================ #
# base class for all genotyping jobs
# ============================================================================================================================ #
# ============================================================================================================================ #
class BaseGenotypingJob(map_reduce.Job):

    def get_output_directory(self):
        c = config.Configuration()
        fastq_file_name = c.fastq_file.split('/')[-1][::-1].split('.')[-1][::-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/genotyping/' + fastq_file_name))

    def get_current_job_directory(self):
        c = config.Configuration()
        # we might be genotyping this sample for various sets of structural variations, keep them separate
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(self.get_output_directory(), self.job_name[:-1], bed_file_name))

    def get_previous_job_directory(self):
        c = config.Configuration()
        # we might be genotyping this sample for various sets of structural variations, keep them separate
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(self.get_output_directory(), self.previous_job_name[:-1], bed_file_name))

# ============================================================================================================================ #
# MapReduce job for exact counting the signature kmers of the sv library in a sample genomew
# ============================================================================================================================ #

class SampleExactKmerCountJob(CountKmersExactJob):

    @staticmethod
    def launch(**kwargs):
        job = SampleExactKmerCountJob(job_name = 'SampleExactKmerCountJob_', previous_job_name = 'MostLikelyBreakPoints_', **kwargs)
        job.execute()

    def check_cli_arguments(self, args):
        # --bed: to specify the set of structural variations we are interested in
        # --fastq: the sample genome we are trying to genotype
        # --threads: the number of processes to fork
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        # load the kmers for this set of structural variations
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as kmers_file:
            self.kmers = {}
            self.tracks = {}
            tracks = json.load(kmers_file)
            for track in tracks:
                self.tracks[track] = {}
                # TODO: change this later
                self.tracks[track]['break_points'] = {}
                self.tracks[track]['break_points'][tracks[track]['most_likey']] = tracks[track]['most_likely']
                for kmer in tracks[track]['most_likely']['novel_kmers']:
                    self.kmers[kmer] = 0
            print('counting signature kmers for', len(tracks), 'tracks totalling', len(self.kmers), 'kmers')
        # cache the break points locally
        with open(os.path.join(self.get_current_job_directory(), 'break_points.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, sort_keys = True, indent = 4)
        # dummy, avoid overrding extra methods
        for index in range(0, self.num_threads):
            self.batch[index] = {}

    def reduce(self):
        c = config.Configuration()
        print('reducing results ...')
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
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4)

    def get_output_directory(self):
        c = config.Configuration()
        fastq_file_name = c.fastq_file.split('/')[-1][::-1].split('.')[-1][::-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/genotyping/' + fastq_file_name))

    def get_current_job_directory(self):
        c = config.Configuration()
        # we might be genotyping this sample for various sets of structural variations, keep them separate
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(self.get_output_directory(), self.job_name[:-1], bed_file_name))

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/' + str(c.ksize) + '/', self.previous_job_name[:-1]))

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job to genotype a genome
# We need to have generated a counttable for the genome we are interested in beforehand
# How to get this far?
# Step 1: get the break points all of which's kmers are found in CHM1, these are the possible candidates for the structural
# variation's boundaries -> BreakPointJob
# Step 2: Find a set of novel kmers for each break point that can be used to indentify it. khmer never underestimates counts so
# if a kmer comes with a count of zero in reference genome, we can be sure that it is really novel -> NovelKmerJob
# Step 3: khmer may report oversetimated counts for these break points so we need to count them exactly again. This is necessary
# for a reliable likelihood model -> CountKmersExactJob
# Step 4: With exact kmer counts available, we can find the most likely break points for each event in our library -> MostLikelyBreakPointsJob
# Step 6: For every given genome we will need to calculate the coverage depth ->
# Step 5: Given a sample genome, try to genotype the structural variations using the likelihood model and signatures gathered
# above -> GenotypingJob (this one)
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenotypingJob(BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = GenotypingJob(job_name = 'Genotyping_', previous_job_name = 'SampleExactKmerCountJob_')
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --coverage option to specify the read depth
        # --std option to specify the coverage standard deviation
        # --fastq to specify the sample being genotyped (this job only needs the name)
        # --bed to indicate the set of structural variations being considered
        pass

    # needs to know the most likely breakpoint and its kmers, MostLikelyBreakPointsJob includes that information
    def load_inputs(self):
        # each event is a structural variation with its most likely breakpoints
        self.tracks = {}
        path = os.path.join(self.get_previous_job_directory(), 'break_points.json')
        with open(path, 'r') as json_file:
            self.tracks = json.load(json_file)
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            self.kmers = json.load(json_file)
        for index in range(0, self.num_threads):
            self.batch[index] = {}
        # round-robin events between available processes
        index = 0
        for track in self.tracks:
            self.batch[index][track] = self.tracks[track]
            index = index + 1
            if index == self.num_threads:
                index = 0

    def transform(self, track, track_name):
        c = config.Configuration()
        likelihood = {}
        distribution = {
            '(1, 1)': statistics.NormalDistribution(mean = c.coverage, std = c.std),
            '(1, 0)': statistics.NormalDistribution(mean = c.coverage / 2, std = c.std),
            '(0, 0)': statistics.ErrorDistribution(p = 1.0 / 1000),
        }
        # all kmers remaining in the track are to be considered part of the signature, no matter the number of breakpoints
        break_points = track['break_points']
        novel_kmers = {}
        for break_point in break_points:
            likelihood[break_point] = {}
            for zyg in distribution:
                likelihood[break_point][zyg] = 0
            novel_kmers[break_point] = {}
            for kmer in break_points[break_point]['novel_kmers']:
                # remember we are genotyping, we need the counts in the sample not in the origin of event
                for zyg in distribution:
                    likelihood[break_point][zyg] += distribution[zyg].log_pmf(self.kmers[kmer])
                novel_kmers[break_point][kmer] = self.kmers[kmer]
        # now find the maximum, for each zygosity find the maximum value, then compare
        for break_point in break_points:
            choice = max(likelihood[break_point].items(), key = operator.itemgetter(1))[0]
            likelihood[break_point] = {'zygosity': choice, 'novel_kmers': novel_kmers[break_point]}
        return likelihood

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    #
    GenotypingJob.launch()
    #SampleExactKmerCountJob.launch()
