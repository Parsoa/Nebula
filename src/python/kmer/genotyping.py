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
    counttable,
    map_reduce,
    statistics,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import plotly.offline as plotly
import plotly.graph_objs as graph_objs

import numpy

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
# MapReduce job for exact counting the signature kmers of the sv library in a sample genome
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
        self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            tracks = json.load(json_file)
        # round-robin events between available processes
        n = 0
        for track in tracks:
            index = n % c.max_threads 
            if not index in self.batch:
                self.batch[index] = {}
            self.batch[index][track] = tracks[track]
            n = n + 1
            self.num_threads = min(n, c.max_threads)

    def transform(self, track, track_name):
        c = config.Configuration()
        likelihood = {}
        distribution = {
            '(1, 1)': statistics.NormalDistribution(mean = c.coverage, std = c.std),
            '(1, 0)': statistics.NormalDistribution(mean = c.coverage / 2, std = c.std),
            '(0, 0)': statistics.NormalDistribution(mean = 0, std = c.std)
        }
        inner_distribution = {
            '(0, 0)': statistics.NormalDistribution(mean = c.coverage, std = c.std),
            '(1, 0)': statistics.NormalDistribution(mean = c.coverage / 2, std = c.std),
            '(1, 1)': statistics.NormalDistribution(mean = 0, std = c.std)
        }
        # all kmers remaining in the track are to be considered part of the signature, no matter the number of breakpoints
        with open(track, 'r') as track_file:
            break_points = json.load(track_file)
        for break_point in break_points:
            likelihood[break_point] = {}
            likelihood[break_point]['novel_kmers'] = {}
            likelihood[break_point]['inner_kmers'] = {}
            likelihood[break_point]['2std'] = 0
            likelihood[break_point]['likelihood'] = {}
            likelihood[break_point]['likelihood']['novel'] = {}
            likelihood[break_point]['likelihood']['inner'] = {}
            for zyg in distribution:
                likelihood[break_point]['likelihood']['inner'][zyg] = 0
                likelihood[break_point]['likelihood']['novel'][zyg] = 0
            for kmer in break_points[break_point]['novel_kmers']:
                count = self.counts_provider.get_kmer_count(kmer)
                if count > c.coverage + 3 * c.std:
                    continue
                likelihood[break_point]['novel_kmers'][kmer] = count
                for zyg in distribution:
                    likelihood[break_point]['likelihood']['novel'][zyg] += distribution[zyg].log_pmf(count)
            for kmer in break_points[break_point]['inner_kmers']:
                count = self.counts_provider.get_kmer_count(kmer)
                if count > c.coverage + 3 * c.std:
                    continue
                likelihood[break_point]['inner_kmers'][kmer] = count
                for zyg in distribution:
                    likelihood[break_point]['likelihood']['inner'][zyg] += inner_distribution[zyg].log_pmf(count)
            for zyg in distribution:
                likelihood[break_point]['likelihood']['inner'][zyg] -= len(likelihood[break_point]['inner_kmers']) * distribution[zyg].log_pmf(distribution[zyg].mean)
                likelihood[break_point]['likelihood']['novel'][zyg] -= len(likelihood[break_point]['novel_kmers']) * distribution[zyg].log_pmf(distribution[zyg].mean)
                likelihood[break_point]['likelihood']['inner'][zyg] = abs(likelihood[break_point]['likelihood']['inner'][zyg])
                likelihood[break_point]['likelihood']['novel'][zyg] = abs(likelihood[break_point]['likelihood']['novel'][zyg])
            inner_choice, inner_cost = min(likelihood[break_point]['likelihood']['inner'].items(), key = operator.itemgetter(1))
            novel_choice, novel_cost = min(likelihood[break_point]['likelihood']['novel'].items(), key = operator.itemgetter(1))
            likelihood[break_point]['genotype'] = {}
            # too few kmers in to do anything useful
            if len(likelihood[break_point]['inner_kmers']) <= 5 and len(likelihood[break_point]['novel_kmers']) <= 5:
                choice = '(2, 2)'
            # enough novel kmers but not enough inner kmers
            elif len(likelihood[break_point]['inner_kmers']) <= 5:
                inner_choice = '(3, 3)'
                choice = novel_choice
            # enough inner kmers but not enough novel kmers
            elif len(likelihood[break_point]['novel_kmers']) <= 5:
                novel_choice = '(4, 4)'
                choice = inner_choice
            else:
                inner_cost = inner_cost / len(likelihood[break_point]['inner_kmers'])
                novel_cost = novel_cost / len(likelihood[break_point]['novel_kmers'])
                choice = inner_choice if inner_cost < novel_cost else novel_choice
            likelihood[break_point]['genotype']['inner'] = inner_choice
            likelihood[break_point]['genotype']['novel'] = novel_choice
            likelihood[break_point]['genotype']['consensus'] = choice
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
        #threshold = []
        #for track in output:
            #with open(output[track], 'r') as json_file:
                #break_points = json.load(json_file)
                #for break_point in break_points:
                    #if 'threshold' in break_points[break_point]:
                        #threshold.append(break_points[break_point]['threshold'])
        #data = [graph_objs.Histogram(x = threshold, xbins = dict(start = 0, size = 5))]
        #plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), 'threshold.html'))
        if c.resume_from_reduce:
            true_positive = bed.read_tracks(os.path.join(self.get_current_job_directory(), 'true_positive.bed'))
            true_negative = bed.read_tracks(os.path.join(self.get_current_job_directory(), 'true_negative.bed'))
            false_positive = bed.read_tracks(os.path.join(self.get_current_job_directory(), 'false_positive.bed'))
            false_negative = bed.read_tracks(os.path.join(self.get_current_job_directory(), 'false_negative.bed'))
            true = {}
            true.update(true_positive)
            true.update(true_negative)
            false = {}
            false.update(false_positive)
            false.update(false_negative)
            b = []
            ratios = []
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            bed_file.write('chrom\tstart\tend\tkmers\tgenotype\t0,0\t1,0\t1,1\tcorrect\n')
            for track in output:
                with open(output[track], 'r') as track_file:
                    break_points = json.load(track_file)
                    chrom = track.split('_')[0]
                    begin = track.split('_')[1]
                    end = track.split('_')[2]
                    for break_point in break_points:
                        inner_kmers = len(break_points[break_point]['inner_kmers'])
                        novel_kmers = len(break_points[break_point]['novel_kmers'])
                        inner_zyg00 = break_points[break_point]['likelihood']['inner']['(0, 0)']
                        inner_zyg10 = break_points[break_point]['likelihood']['inner']['(1, 0)']
                        inner_zyg11 = break_points[break_point]['likelihood']['inner']['(1, 1)']
                        novel_zyg00 = break_points[break_point]['likelihood']['novel']['(0, 0)']
                        novel_zyg10 = break_points[break_point]['likelihood']['novel']['(1, 0)']
                        novel_zyg11 = break_points[break_point]['likelihood']['novel']['(1, 1)']
                        inner = break_points[break_point]['genotype']['inner']
                        novel = break_points[break_point]['genotype']['novel']
                        consensus = break_points[break_point]['genotype']['consensus']
                        bed_file.write(chrom + '\t' + begin + '\t' + end + '\t' + str(novel_kmers) + '\t' + str(inner_kmers) + '\t' +
                                str(inner_zyg00 + novel_zyg00) + '\t' +
                                str(inner_zyg10 + novel_zyg10) + '\t' +
                                str(inner_zyg11 + novel_zyg11) + '\t' +
                                str(consensus) + '\t' +
                                '\n')
        if c.resume_from_reduce:
            p = numpy.corrcoef(b, ratios)
            print(p)

    def plot(self, b):
        zyg11_t = list(map(lambda x: x[1], list(filter(lambda l: l[0] == '(1, 1)' and l[2], b))))
        zyg11_f = list(map(lambda x: x[1], list(filter(lambda l: l[0] == '(1, 1)' and not l[2], b))))
        zyg10_t = list(map(lambda x: x[1], list(filter(lambda l: l[0] == '(1, 0)' and l[2], b))))
        zyg10_f = list(map(lambda x: x[1], list(filter(lambda l: l[0] == '(1, 0)' and not l[2], b))))
        zyg00_t = list(map(lambda x: x[1], list(filter(lambda l: l[0] == '(0, 0)' and l[2], b))))
        zyg00_f = list(map(lambda x: x[1], list(filter(lambda l: l[0] == '(0, 0)' and not l[2], b))))
        trace0_t = graph_objs.Scatter(x = zyg11_t, y = [1] * len(zyg11_t), mode = 'markers', marker = dict(color = 'rgb(66, 134, 244)'))
        trace0_f = graph_objs.Scatter(x = zyg11_f, y = [1] * len(zyg11_f), mode = 'markers', marker = dict(color = 'rgb(244, 65, 65)'))
        trace1_t = graph_objs.Scatter(x = zyg10_t, y = [2] * len(zyg10_t), mode = 'markers', marker = dict(color = 'rgb(66, 134, 244)'))
        trace1_f = graph_objs.Scatter(x = zyg10_f, y = [2] * len(zyg10_f), mode = 'markers', marker = dict(color = 'rgb(244, 65, 65)'))
        trace2_t = graph_objs.Scatter(x = zyg00_t, y = [3] * len(zyg00_t), mode = 'markers', marker = dict(color = 'rgb(66, 134, 244)'))
        trace2_f = graph_objs.Scatter(x = zyg00_f, y = [3] * len(zyg00_f), mode = 'markers', marker = dict(color = 'rgb(244, 65, 65)'))
        data = [trace0_t, trace0_f, trace1_t, trace1_f, trace2_t, trace2_f]
        path = os.path.join(self.get_current_job_directory(), 'scatter')
        plotly.plot(data, filename = path)

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/' + str(c.ksize) + '/', self.previous_job_name[:-1]))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenotypingAnalysisJob(BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = GenotypingAnalysisJob(job_name = 'Genotyping_', previous_job_name = 'Genotyping_', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = 1

    def load_inputs(self):
        self.true_positive = {}
        self.true_negative = {}
        self.false_positive = {}
        self.false_negative = {}
        classes = ['true_positive', 'true_negative', 'false_positive', 'false_negative']
        for c in classes:
            setattr(self, c, {})
            with open(os.path.join(self.get_current_job_directory(), c + '.bed'), 'r') as bed_file:
                for track in bed.parse_bed_file(bed_file):
                    track_name = track[0] + '_' + track[1] + '_' + track[2]
                    getattr(self, c)[track_name] = track
        for c in classes:
            self.plot_kmer_count_variation(c, getattr(self, c))
        exit()

    def plot_kmer_count_variation(self, batch_name, batch):
        variance = {'inner': [], 'novel': []}
        for track_name in batch:
            path = os.path.join(self.get_current_job_directory(), 'genotype_' + track_name + '.json')
            if not os.path.exists(path):
                print(red('not found:', path))
                continue
            with open(path, 'r') as json_file:
                break_points = json.load(json_file)
                for break_point in break_points:
                    inner_kmers = list(map(lambda x: x[1], break_points[break_point]['inner_kmers'].items()))
                    novel_kmers = list(map(lambda x: x[1], break_points[break_point]['novel_kmers'].items()))
                    if len(inner_kmers):
                        variance['inner'].append(statistics.variance(inner_kmers))
                    if len(novel_kmers):
                        variance['novel'].append(statistics.variance(novel_kmers))
        v = statistics.variance(variance['inner'])
        m = statistics.mean(variance['inner'])
        data = [graph_objs.Histogram(x = variance['inner'], xbins = dict(start = 0, size = 5, end = 1000), text = 'Mean:' + str(m) + ' STD: ' + str(math.sqrt(v)))]
        plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), batch_name + '_inner.html'), auto_open = False)
        v = statistics.variance(variance['novel'])
        m = statistics.mean(variance['novel'])
        data = [graph_objs.Histogram(x = variance['novel'], xbins = dict(start = 5, size = 5, end = 1000), text = 'Mean:' + str(m) + ' STD: ' + str(math.sqrt(v)))]
        plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), batch_name + '_novel.html'), auto_open = False)

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    if c.job == 'GenotypingJob':
        GenotypingJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'GenotypingAnalysisJob':
        GenotypingAnalysisJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'SampleExactKmerCountingJob':
        SampleExactKmerCountJob.launch(resume_from_reduce = c.resume_from_reduce)
