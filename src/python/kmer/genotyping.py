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

from kmer import (
    bed,
    config,
    map_reduce,
    statistics,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

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

class SampleExactKmerCountingJob(map_reduce.BaseExactCountingJob):

    @staticmethod
    def launch(**kwargs):
        job = SampleExactKmerCountJob(job_name = 'SampleExactKmerCountingJob_', previous_job_name = 'MostLikelyBreakPoints_', **kwargs)
        job.execute()

    def check_cli_arguments(self, args):
        # --bed: to specify the set of structural variations we are interested in
        # --fastq: the sample genome we are trying to genotype
        # --threads: the number of processes to fork
        pass

    def load_inputs(self):
        # load the kmers for this set of structural variations
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        self.kmers = {}
        self.tracks = {}
        with open(path, 'r') as tracks_file:
            paths = json.load(tracks_file)
            for track in paths:
                self.tracks[track] = {}
                with open(paths[track], 'r') as json_file:
                    break_points = json.load(json_file)
                    self.tracks[track]['break_points'] = break_points
                    for break_point in break_points:
                        for kmer in break_points[break_point]['novel_kmers']:
                            self.kmers[kmer] = 0
        print('counting signature kmers for', len(self.tracks), 'tracks totalling', len(self.kmers), 'kmers')
        # dummy, avoid overrding extra methods
        for index in range(0, self.num_threads):
            self.batch[index] = {}

    def reduce(self):
        c = config.Configuration()
        kmers = self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4)
        # update break_points
        output = {}
        for track in self.tracks:
            for break_point in self.tracks[track]['break_points']:
                for kmer in self.tracks[track]['break_points'][break_point]['novel_kmers']:
                    break_points[break_point]['novel_kmers'][kmer].pop('count', None)
                    break_points[break_point]['novel_kmers'][kmer].pop('actual_count', None)
                    break_points[break_point]['novel_kmers'][kmer]['sample_count'] = kmers[kmer]
            path = os.path.join(self.get_current_job_directory(), 'sample_exact_counts_' + track + '.json')
            output[track] = path
            with open(path, 'w') as json_file:
                json.dump(self.tracks[track], json_file, sort_keys = True, indent = 4)
        # dump the kmer counts
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)

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
    def launch(**kwargs):
        job = GenotypingJob(job_name = 'Genotyping_', previous_job_name = 'MostLikelyBreakPointsJob_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def check_cli_arguments(self, args):
        # --coverage option to specify the read depth
        # --std option to specify the coverage standard deviation
        # --fastq to specify the sample being genotyped (this job only needs the name)
        # --bed to indicate the set of structural variations being considered (only needs the name)
        pass

    # needs to know the most likely breakpoint and its kmers, MostLikelyBreakPointsJob includes that information
    def load_inputs(self):
        # each event is a structural variation with its most likely breakpoints
        self.tracks = {}
        self.counts_provider = None
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            self.tracks = json.load(json_file)
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
            #'(0, 0)': statistics.ErrorDistribution(p = 1.0 / 100),
            '(0, 0)': statistics.NormalDistribution(mean = 0, std = 5)
        }
        # all kmers remaining in the track are to be considered part of the signature, no matter the number of breakpoints
        with open(track, 'r') as track_file:
            print(track)
            break_points = json.load(track_file)
        novel_kmers = {}
        for break_point in break_points:
            likelihood[break_point] = {}
            likelihood[break_point]['kmers'] = {}
            likelihood[break_point]['2std'] = 0
            likelihood[break_point]['zyg'] = {}
            for zyg in distribution:
                likelihood[break_point]['zyg'][zyg] = 0
            for kmer in break_points[break_point]['novel_kmers']:
                count = self.get_kmer_count(kmer, self.index, False)
                likelihood[break_point]['kmers'][kmer] = count
                # remember we are genotyping, we need the counts in the sample not in the origin of event
                for zyg in distribution:
                    likelihood[break_point]['zyg'][zyg] += distribution[zyg].log_pmf(count)

                if abs(count - c.coverage) < 2 * c.std:
                    likelihood[break_point]['2std'] += 1
        # now find the maximum, for each zygosity find the maximum value, then compare
        for break_point in break_points:
            choice = max(likelihood[break_point]['zyg'].items(), key = operator.itemgetter(1))[0]
            likelihood[break_point]['genotype'] = choice#, 'novel_kmers': novel_kmers[break_point]}
        path = os.path.join(self.get_current_job_directory(), 'genotype_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(likelihood, json_file, sort_keys = True, indent = 4)
        return path

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
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            for track in output:
                with open(output[track], 'r') as track_file:
                    break_points = json.load(track_file)
                    chrom = track.split('_')[0]
                    begin = track.split('_')[1]
                    end = track.split('_')[2]
                    for break_point in break_points:
                        std2 = break_points[break_point]['2std']
                        kmers = len(break_points[break_point]['kmers'])
                        zyg00 = break_points[break_point]['zyg']['(0, 0)']
                        zyg10 = break_points[break_point]['zyg']['(1, 0)']
                        zyg11 = break_points[break_point]['zyg']['(1, 1)']
                        genotype = break_points[break_point]['genotype']
                        bed_file.write(chrom + '\t' + begin + '\t' + end + '\t' + str(kmers) + '\t' + str(std2) + '\t' + genotype + '\t' + str(zyg00) + '\t' + str(zyg10) + '\t' + str(zyg11) + '\n')

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/' + str(c.ksize) + '/', self.previous_job_name[:-1]))

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    if c.job == 'GenotypingJob':
        GenotypingJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'SampleExactKmerCountingJob':
        SampleExactKmerCountJob.launch(resume_from_reduce = c.resume_from_reduce)
