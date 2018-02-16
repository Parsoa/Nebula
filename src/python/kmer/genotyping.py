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

class SampleExactKmerCountJob(map_reduce.BaseExactCountingJob):

    @staticmethod
    def launch(**kwargs):
        job = SampleExactKmerCountJob(job_name = 'SampleExactKmerCountJob_', previous_job_name = 'MostLikelyBreakPoints_', **kwargs)
        job.execute()

    def check_cli_arguments(self, args):
        # --bed: to specify the set of structural variations we are interested in
        # --fastq: the sample genome we are trying to genotype
        # --threads: the number of processes to fork
        pass

    def load_inputs(self):
        # load the kmers for this set of structural variations
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        # we are only interested in the most likely break points from each event
        with open(path, 'r') as kmers_file:
            self.kmers = {}
            self.tracks = {}
            tracks = json.load(kmers_file)
            for track in tracks:
                self.tracks[track] = {}
                self.tracks[track]['break_points'] = tracks[track]['most_likely']
                for break_point in tracks[track]['most_likely']:
                    for kmer in tracks[track]['most_likely'][break_point]['novel_kmers']:
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
        kmers = self.merge_counts()
        # update break_points
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as kmers_file:
            tracks = json.load(kmers_file)
            for track in tracks:
                tracks[track]['break_points'] = tracks[track]['most_likely']
                tracks[track].pop('most_likely', None)
        with open(os.path.join(self.get_current_job_directory(), 'break_points.json'), 'r') as json_file:
            tracks = json.load(json_file)
        for track in tracks:
            for break_point in tracks[track]['break_points']:
                for kmer in tracks[track]['break_points'][break_point]['novel_kmers']:
                     tracks[track]['break_points'][break_point]['novel_kmers'][kmer] = kmers[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'break_points.json'), 'w') as json_file:
            json.dump(tracks, json_file, sort_keys = True, indent = 4)
        # dump the kmer counts
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
        # --bed to indicate the set of structural variations being considered (only needs the name)
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
            '(0, 0)': statistics.ErrorDistribution(p = 1.0 / 100),
        }
        # all kmers remaining in the track are to be considered part of the signature, no matter the number of breakpoints
        break_points = track['break_points']
        novel_kmers = {}
        for break_point in break_points:
            likelihood[break_point] = {}
            likelihood[break_point]['kmers'] = len(break_points[break_point]['novel_kmers'])
            likelihood[break_point]['2std'] = 0
            likelihood[break_point]['zyg'] = {}
            for zyg in distribution:
                likelihood[break_point]['zyg'][zyg] = 0
            for kmer in break_points[break_point]['novel_kmers']:
                # remember we are genotyping, we need the counts in the sample not in the origin of event
                for zyg in distribution:
                    likelihood[break_point]['zyg'][zyg] += distribution[zyg].log_pmf(self.kmers[kmer])

                if abs(self.kmers[kmer] - c.coverage) < 2 * c.std:
                    likelihood[break_point]['2std'] += 1
        # now find the maximum, for each zygosity find the maximum value, then compare
        for break_point in break_points:
            choice = max(likelihood[break_point]['zyg'].items(), key = operator.itemgetter(1))[0]
            likelihood[break_point]['max'] = choice#, 'novel_kmers': novel_kmers[break_point]}
        return likelihood

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
                chrom = track.split('_')[0]
                begin = track.split('_')[1]
                end = track.split('_')[2]
                for break_point in output[track]:
                    kmers = output[track][break_point]['kmers']
                    std2 = output[track][break_point]['2std']
                    choice = output[track][break_point]['max']
                    zyg00 = output[track][break_point]['zyg']['(0, 0)']
                    zyg10 = output[track][break_point]['zyg']['(1, 0)']
                    zyg11 = output[track][break_point]['zyg']['(1, 1)']
                    bed_file.write(chrom + '\t' + begin + '\t' + end + '\t' + str(kmers) + '\t' + str(std2) + '\t' + choice + '\t' + str(zyg00) + '\t' + str(zyg10) + '\t' + str(zyg11) + '\n')

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    #
    GenotypingJob.launch()
    #SampleExactKmerCountJob.launch(resume_from_reduce = True)
