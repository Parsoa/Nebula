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
    training,
    counttable,
    map_reduce,
    statistics,
)

from kmer.sv import StructuralVariation, Inversion, Deletion, SNP
from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import pybedtools

import tensorflow as tf
import plotly.offline as plotly
import plotly.graph_objs as graph_objs

# ============================================================================================================================ #
# ============================================================================================================================ #
# base class for all genotyping jobs
# ============================================================================================================================ #
# ============================================================================================================================ #

class BaseTensorflowJob(training.BaseTrainingJob):

    def get_output_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../training/' + bed_file_name + '/Tensorflow/'))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class TensorflowTrainingJob(BaseTensorflowJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = TensorflowTrainingJob(job_name = 'TensorflowTrainingJob_', previous_job_name = '', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform_class_labels(self, Y):
        return list(map(lambda y: [1.0, 0.0, 0.0] if y == 0 else [0.0, 1.0, 0.0] if y == 1 else [0.0, 0.0, 1.0], Y))
        #return map(list(lambda y: [y], Y))

    def transform(self, track, track_name):
        c = config.Configuration()
        path = os.path.join(self.get_current_job_directory(), track_name)
        if not os.path.exists(path):
            os.makedirs(path)
        tf.reset_default_graph()
        learning_rate = 0.1
        inner_kmers, novel_kmers, _1, _2 = self.extract_kmers(track, track_name)
        X, Y = self.simulate_kmer_counts(track, track_name, 20)
        num_features = len(inner_kmers) + len(novel_kmers)
        x = tf.placeholder(tf.float32, [None, num_features], name = 'x')
        y = tf.placeholder(tf.float32, [None, 3], name = 'y')
        Y = self.transform_class_labels(Y)
        print(len(X), len(Y))
        json_print({'y': Y})
        model = self.create_network(num_features, 3, 1, x)
        #model = tf.Print(model, [model])
        loss = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits = model, labels = y))
        optimizer = tf.train.AdamOptimizer(learning_rate = learning_rate).minimize(loss)
        init = tf.global_variables_initializer()
        saver = tf.train.Saver()
        with tf.Session() as session:
            session.run(init)
            # how many?
            for epoch in range(10):
                o, c = session.run([optimizer, loss], feed_dict = {x: X, y: Y})
            r = model.eval({x: X, y: Y})
            print(r)
            #correct_prediction = tf.equal(tf.argmax(model, 1), tf.argmax(Y, 1))
            # Calculate accuracy
            #accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
            #print(model)
            print(green("Accuracy:", accuracy.eval({x: X, y: Y})))
            #model_path = os.path.join(path, 'model.tf')
            #save_path = saver.save(session, model_path)
            print(save_path)
        print(cyan(track_name), green(num_features, 'kmer'), blue(len(X[0]), 'features'))
        with open(os.path.join(path, 'kmers.json'), 'w') as json_file:
            json.dump({'inner_kmers': inner_kmers, 'novel_kmers': novel_kmers, 'features': num_features}, json_file, sort_keys = True, indent = 4)
        return path
        tf.reset_default_graph()
        with tf.Session() as session:
            saver = tf.train.import_meta_graph(save_path + '.index')
            saver.restore(session, tf.train.latest_checkpoint(os.path.join(self.get_current_job_directory(), track_name)))
            x = tf.get_default_graph().get_tensor_by_name('x:0')
            print(x)
            #print(cyan(track_name), green(num_features, 'kmer'), blue(len(X[0]), 'features'), purple(x))
        return path

    def create_network(self, num_features, num_classes, num_layers, x):
        weights = []
        for i in range(num_layers):
            weights.append(tf.Variable(tf.random_normal([num_features, num_features])))
        weights.append(tf.Variable(tf.random_normal([num_features, num_classes])))
        biases = []
        for i in range(num_layers):
            biases.append(tf.Variable(tf.random_normal([num_features])))
        biases.append(tf.Variable(tf.random_normal([num_classes])))
        layers = []
        layers.append(tf.add(tf.matmul(x, weights[0]), biases[0]))
        layers[0] = tf.nn.sigmoid(layers[0])
        for i in range(1, num_layers):
            layers.append(tf.add(tf.matmul(layers[i - 1], weights[i]), biases[i]))
            layers[i] = tf.sigmoid.relu(layers[i])
        # Output layer with linear activation
        out_layer = tf.add(tf.matmul(layers[num_layers - 1], weights[num_layers]), biases[num_layers], name = 'output')
        print(out_layer)
        return out_layer

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

class TensorflowGenotypingJob(BaseTensorflowJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = TensorflowGenotypingJob(job_name = 'TensorflowGenotypingJob_', previous_job_name = 'TensorflowTrainingJob_', **kwargs)
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
        init = tf.global_variables_initializer()
        graph = tf.Graph()
        session = tf.Session(graph = graph)
        with graph.as_default():
            print(cyan(track))
            saver = tf.train.import_meta_graph(os.path.join(track, 'model.tf.meta'))
            saver.restore(session, tf.train.latest_checkpoint(track))
            inner_kmers, novel_kmers, num_features = self.count_kmers_in_sample(track)
            num_classes = 3
            num_features = len(inner_kmers) + len(novel_kmers)
            x = graph.get_tensor_by_name('x:0')
            features = []
            features += list(map(lambda t: t[1] / (1.0 * c.coverage), sorted(list(map(lambda kmer: (kmer, inner_kmers[kmer]), inner_kmers)), key = operator.itemgetter(0))))
            features += list(map(lambda t: t[1] / (1.0 * c.coverage), sorted(list(map(lambda kmer: (kmer, novel_kmers[kmer]), novel_kmers)), key = operator.itemgetter(0))))
            print(green(len(features), 'kmers'), blue(num_features, 'features'), cyan(x))
            try:
                p = session.run(graph.get_tensor_by_name('output:0'), feed_dict = {x: [features]})
            except Exception as e:
                print(red('Tensor shape mismatch'))
                return None
            choice = tf.argmax(p, 1)
            genotype = '(0, 0)' if choice == 0 else '(1, 0)' if choice == 1 else '(1, 1)'
            path = os.path.join(self.get_current_job_directory(), 'genotype_' + track_name + '.json')
            with open(path, 'w') as json_file:
                json.dump({'genotype': genotype, 'inner_kmers': inner_kmers, 'novel_kmers': novel_kmers}, json_file, indent = 4) 
            return path

    def count_kmers_in_sample(self, path):
        with open(os.path.join(path, 'kmers.json'), 'r') as json_file:
            k = json.load(json_file)
            inner_kmers = k['inner_kmers']
            novel_kmers = k['novel_kmers']
            for kmer in novel_kmers:
                count = self.counts_provider.get_kmer_count(kmer)
                novel_kmers[kmer] = count / float(c.coverage)
            for kmer in inner_kmers:
                count = self.counts_provider.get_kmer_count(kmer)
                inner_kmers[kmer] = count / float(c.coverage)
            return inner_kmers, novel_kmers, k['features']

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
                with open(output[track], 'r') as json_file:
                    print(track)
                    payload = json.load(json_file) 
                    chrom = track.split('_')[0]
                    begin = track.split('_')[1]
                    end = track.split('_')[2]
                    inner_kmers = len(payload['inner_kmers'])
                    novel_kmers = len(payload['novel_kmers'])
                    genotype = payload['genotype']
                    bed_file.write(chrom + '\t' + begin + '\t' + end + '\t' + str(novel_kmers) + '\t' + str(inner_kmers) + '\t' +
                            genotype + '\n')

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
    if c.job == 'TensorflowTrainingJob':
        TensorflowTrainingJob.launch(resume_from_reduce = c.resume_from_reduce)
    if c.job == 'TensorflowGenotypingJob':
        TensorflowGenotypingJob.launch(resume_from_reduce = c.resume_from_reduce)
