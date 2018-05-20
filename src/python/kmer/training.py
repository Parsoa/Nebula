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
import cPickle
import argparse
import operator
import traceback
import subprocess

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

from sklearn import svm
from sklearn import tree

# ============================================================================================================================ #
# ============================================================================================================================ #
# base class for all genotyping jobs
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseDecisionTreeJob(map_reduce.Job):

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../training/' + bed_file_name + '/'))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class DecisionTreeTrainingJob(BaseDecisionTreeJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = DecisionTreeTrainingJob(job_name = 'DecisionTreeTrainingJob_', previous_job_name = '', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'merge.json'), 'r') as json_file:
            tracks = json.load(json_file)
        n = 0
        for track in tracks:
            tokens = track.split('_')
            start = int(tokens[1])
            end = int(tokens[2])
            chrom = tokens[0]
            #if chrom != 'chr22':
                #continue
            #if start != 20950622:
                #continue
            #if end != 20953839:
                #continue
            index = n % c.max_threads 
            if not index in self.batch:
                self.batch[index] = {}
            self.batch[index][track] = tracks[track]
            print(blue('assigned', track, 'to', index))
            n = n + 1
            self.num_threads = min(n, c.max_threads)

    def transform(self, track, track_name):
        start = time.time()
        c = config.Configuration()
        # load kmers and chosen breakpoint
        self.inner_kmers = {}
        self.novel_kmers = {}
        with open(track) as json_file:
            break_points = json.load(json_file)
            for break_point in break_points:
                for kmer in break_points[break_point]['inner_kmers']:
                    canon = get_canonical_kmer_representation(kmer)
                    self.inner_kmers[canon] = 0
                for kmer in break_points[break_point]['novel_kmers']:
                    canon = get_canonical_kmer_representation(kmer)
                    self.novel_kmers[canon] = 0
                start_offset = int(break_point.replace('(', '').replace(')', '').split(',')[0])
                end_offset = int(break_point.replace('(', '').replace(')', '').split(',')[1])
        # extract sequence and apply event
        tokens = track_name.split('_')
        start = int(tokens[1]) + start_offset
        end = int(tokens[2]) + end_offset
        chrom = extract_chromosome(tokens[0])
        offset = 1000
        base = chrom[start - offset - 1 : end + offset - 1]
        event = base[:offset] + base[-offset:]
        random.seed(c.seed)
        X = []
        Y = []
        path = os.path.join(self.get_current_job_directory(), track_name)
        if len(self.novel_kmers) + len(self.inner_kmers) < 5:
            print(red('skipping'), track_name)
            return None
        print(track_name, len(self.novel_kmers) + len(self.inner_kmers), 'kmers')
        if not os.path.exists(path):
            os.makedirs(path)
        path = os.path.join(self.get_current_job_directory(), track_name)
        with open(os.path.join(path, 'kmers.json'), 'w') as json_file:
            json.dump({'inner_kmers': self.inner_kmers, 'novel_kmers': self.novel_kmers}, json_file, sort_keys = True, indent = 4)
        for i in range(0, 20):
            self.kmers = {}
            self.kmers.update({'inner': {get_canonical_kmer_representation(kmer): 0 for kmer in self.inner_kmers}})
            self.kmers.update({'novel': {get_canonical_kmer_representation(kmer): 0 for kmer in self.novel_kmers}})
            coverage = random.randint(7, 50)
            zygosity = random.randint(0, 2)
            #print('iteration', blue(i), 'coverage', green(coverage), 'zygosity', cyan(zygosity))
            if zygosity == 0:
                strands = [base, base]
            elif zygosity == 1:
                strands = [base, event]
            else:
                strands = [event, event]
            for s in range(0, 2):
                base_name = track_name + '_' + str(coverage) + 'x_' + str(i) + '.' + str(s)
                base_path = os.path.join(self.get_current_job_directory(), track_name, base_name)
                self.export_fasta(strands[s], base_path, track_name)
                self.export_fastq(strands[s], base_path, track_name, coverage)
                self.count_kmers_in_fastq(base_path)
                # cleanup
                os.remove(base_path + '.fa')
                os.remove(base_path + '.0.fq')
                os.remove(base_path + '.1.fq')
            # normalize kmer counts
            #print(blue('======================================================================'))
            #print(zygosity, coverage)
            #json_print(self.kmers)
            # update training date
            features = []
            features += list(map(lambda t: t[1] / (2.0 * coverage), sorted(list(map(lambda kmer: (kmer, self.kmers['inner'][kmer]), self.kmers['inner'])), key = operator.itemgetter(0))))
            features += list(map(lambda t: t[1] / (2.0 * coverage), sorted(list(map(lambda kmer: (kmer, self.kmers['novel'][kmer]), self.kmers['novel'])), key = operator.itemgetter(0))))
            X.append(features)
            #print(len(features), 'festures')
            Y.append(zygosity)
        #clf = tree.DecisionTreeClassifier()
        clf = svm.SVC()
        clf.fit(X, Y)
        #dot_data = tree.export_graphviz(clf, out_file = os.path.join(path, 'tree.dot'))
        with open(os.path.join(path, 'decision_tree.pkl'), 'wb') as pickle_file:
            cPickle.dump(clf, pickle_file)
        #rules, n = self.export_decision_tree_rules(clf)
        #with open(os.path.join(path, 'decision_tree.txt'), 'w') as rules_file:
            #rules_file.write(rules)
        # export kmers
        print('took', blue(time.time() - start))
        return path

    def export_fasta(self, sequence, path, track):
        c = config.Configuration()
        with open(path + '.fa', 'w') as fasta_file:
            w = 50
            l = len(sequence)
            num_lines = l / w
            fasta_file.write('>' + track + '\n')
            for i in range(0, num_lines):
                fasta_file.write(sequence[i * w : (i + 1) * w] + '\n')
            if l % w != 0:
                r = w - l % w
                s = 'N' * r
                fasta_file.write(sequence[num_lines * w :] + s + '\n')

    def export_fastq(self, sequence, path, track, coverage):
        FNULL = open(os.devnull, 'w')
        c = config.Configuration()
        num_reads = len(sequence) * coverage / 100
        fasta = os.path.join(path + '.fa')
        fastq_0 = os.path.join(path + '.0.fq')
        fastq_1 = os.path.join(path + '.1.fq')
        command = "wgsim -d400 -N{} -1100 -2100 {} {} {}".format(num_reads, fasta, fastq_0, fastq_1)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def count_kmers_in_fastq(self, path):
        c = config.Configuration()
        for i in range(0, 2):
            for seq, name in parse_fastq(path + '.' + str(i) + '.fq'):
                kmers = extract_canonical_kmers(31, seq)
                for kmer in kmers:
                    if kmer in self.kmers['inner']:
                        self.kmers['inner'][kmer] += kmers[kmer]
                    if kmer in self.kmers['novel']:
                        self.kmers['novel'][kmer] += kmers[kmer]

    def export_decision_tree_rules(self, clf):
        tree_ = clf.tree_
        feature_names = []
        feature_names += list(map(lambda t: t[0], sorted(list(map(lambda kmer: (kmer, self.kmers['inner'][kmer]), self.kmers['inner'])), key = operator.itemgetter(0))))
        feature_names += list(map(lambda t: t[0], sorted(list(map(lambda kmer: (kmer, self.kmers['novel'][kmer]), self.kmers['novel'])), key = operator.itemgetter(0))))
        feature_names = [
            feature_names[i] if i != tree._tree.TREE_UNDEFINED else "undefined!"
            for i in tree_.feature
        ]
        def recurse(node, depth, rules):
            print('depth:', depth)
            indent = "  " * depth
            if tree_.feature[node] != tree._tree.TREE_UNDEFINED:
                name = feature_names[node]
                threshold = tree_.threshold[node]
                rules += "{}if {} <= {}:\n".format(indent, name, threshold)
                rules, r = recurse(tree_.children_left[node], depth + 1, rules)
                rules += "{}else:  # if {} > {}\n".format(indent, name, threshold)
                rules, l = recurse(tree_.children_right[node], depth + 1, rules)
                return rules, r + l
            else:
                rules += "{}return {}\n".format(indent, tree_.value[node])
                return rules, 1
        rules, m = recurse(0, 1, '')
        print('using', green(len(clf.feature_importances_)), 'features')
        return rules, m

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/31/MostLikelyBreakPointsJob/'))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class DecisionTreeGenotypingJob(BaseDecisionTreeJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = DecisionTreeGenotypingJob(job_name = 'DecisionTreeGenotypingJob_', previous_job_name = 'DecisionTreeTrainingJob_', **kwargs)
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

    def load_inputs(self):
        c = config.Configuration()
        self.tracks = {}
        self.counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        path = os.path.join(self.get_previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            self.tracks = json.load(json_file)
        for index in range(0, self.num_threads):
            self.batch[index] = {}
        # round-robin events between available processes
        index = 0
        n = 0
        for track in self.tracks:
            self.batch[index][track] = self.tracks[track]
            index = index + 1
            if index == self.num_threads:
                index = 0
            n += 1

    def transform(self, track, track_name):
        c = config.Configuration()
        with open(os.path.join(track, 'kmers.json'), 'r') as json_file:
            k = json.load(json_file)
            inner_kmers = k['inner_kmers']
            novel_kmers = k['novel_kmers']
            for kmer in novel_kmers:
                count = self.counts_provider.get_kmer_count(kmer)
                novel_kmers[kmer] = count / float(c.coverage)
            for kmer in inner_kmers:
                count = self.counts_provider.get_kmer_count(kmer)
                inner_kmers[kmer] = count / float(c.coverage)
        kmers = {}
        kmers.update({'inner': {get_canonical_kmer_representation(kmer): inner_kmers[kmer] for kmer in inner_kmers}})
        kmers.update({'novel': {get_canonical_kmer_representation(kmer): novel_kmers[kmer] for kmer in novel_kmers}})
        features = []
        features += list(map(lambda t: t[1], sorted(list(map(lambda kmer: (kmer, kmers['inner'][kmer]), kmers['inner'])), key = operator.itemgetter(0))))
        features += list(map(lambda t: t[1], sorted(list(map(lambda kmer: (kmer, kmers['novel'][kmer]), kmers['novel'])), key = operator.itemgetter(0))))
        with open(os.path.join(track, 'decision_tree.pkl'), 'r') as tree_file:
            decision_tree = cPickle.load(tree_file)
        distance = decision_tree.decision_function([features])
        print(distance)
        prediction = decision_tree.predict([features])
        print(prediction)
        genotype = '(0, 0)' if prediction == 0 else '(1, 0)' if prediction == 1 else '(1, 1)'
        path = os.path.join(self.get_current_job_directory(), 'genotype_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({'genotype': genotype, 'inner_kmers': inner_kmers, 'novel_kmers': novel_kmers, 'distance': {'(0, 0)': distance[0][0], '(1, 0)': distance[0][1], '(1, 1)': distance[0][2]}}, json_file, sort_keys = True, indent = 4)
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
        d_00 = []
        d_10 = []
        d_11 = []
        with open(os.path.join(self.get_current_job_directory(), 'merge.bed'), 'w') as bed_file:
            for track in output:
                with open(output[track], 'r') as json_file:
                    print(track)
                    payload = json.load(json_file) 
                    chrom = track.split('_')[0]
                    begin = track.split('_')[1]
                    end = track.split('_')[2]
                    inner_kmers = len(payload['inner_kmers'])
                    novel_kmers = len(payload['novel_kmers'])
                    genotype = payload['genotype']
                    distance_00 = payload['distance']['(0, 0)']
                    if genotype == '(0, 0)':
                        d_00.append(distance_00)
                    distance_10 = payload['distance']['(1, 0)']
                    if genotype == '(1, 0)':
                        d_10.append(distance_10)
                    distance_11 = payload['distance']['(1, 1)']
                    if genotype == '(1, 1)':
                        d_11.append(distance_11)
                    bed_file.write(chrom + '\t' + begin + '\t' + end + '\t' + str(novel_kmers) + '\t' + str(inner_kmers) + '\t' +
                            str(distance_00) + '\t' + str(distance_10) + '\t' + str(distance_11) + '\t' +
                            genotype + '\n')
        data = [graph_objs.Scatter(y = d_00, x = list(range(0, len(d_00))))]
        plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), '00.html'), auto_open = False)
        data = [graph_objs.Scatter(y = d_10, x = list(range(0, len(d_10))))]
        plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), '10.html'), auto_open = False)
        data = [graph_objs.Scatter(y = d_11, x = list(range(0, len(d_11))))]
        plotly.plot(data, filename = os.path.join(self.get_current_job_directory(), '11.html'), auto_open = False)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    if c.job == 'DecisionTreeTrainingJob':
        DecisionTreeTrainingJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'DecisionTreeGenotypingJob':
        DecisionTreeGenotypingJob.launch(resume_from_reduce = c.resume_from_reduce)
