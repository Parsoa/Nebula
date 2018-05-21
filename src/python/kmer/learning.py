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
# ============================================================================================================================ #

class TensorflowTrainingJob(training.BaseTrainingJob):

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
        return map(list(lambda y: [1.0, 0.0, 0.0] if y == 0 else [0.0, 1.0, 0.0] if y == 1 else [0.0, 0.0, 1.0], Y))
        #return map(list(lambda y: [y], Y))

    def transform(self, track, track_name):
        c = config.Configuration()
        path = os.path.join(self.get_current_job_directory(), track_name)
        if not os.path.exists(path):
            os.makedirs(path)
        learning_rate = 0.1
        inner_kmers, novel_kmers = self.extract_kmers(track, track_name)
        X, Y = self.simulate_kmer_counts(track, track_name, 20)
        x = tf.placeholder(tf.float32, [None, num_inputs])
        y = tf.placeholder(tf.float32, [None, 3])
        Y = self.transform_class_labels(Y)
        model = self.create_network(len(inner_kmers) + len(novel_kmers), 3, x)
        cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits = model, labels = y))
        optimizer = tf.train.AdamOptimizer(learning_rate = learning_rate).minimize(cost)
        init = tf.global_variables_initializer()
        saver = tf.train.Saver()
        with tf.Session() as session:
            session.run(init)
            # how many?
            for epoch in range(5):
                o, c = sessions.run([optimizer, cost], feed_dict = {x: X, y: Y})
            correct_prediction = tf.equal(tf.argmax(model, 1), tf.argmax(Y, 1))
            # Calculate accuracy
            #accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
            #print("Accuracy:", accuracy.eval({x: , y: }))
        model_path = os.path.join(path, 'model.tf')
        return path

    def create_network(self, num_inputs, num_classes, x):
        weights = {
            'h1': tf.Variable(tf.random_normal([num_inputs, num_inputs])),
            'h2': tf.Variable(tf.random_normal([num_inputs, num_inputs])),
            'out': tf.Variable(tf.random_normal([num_inputs, num_classes]))
        }
        biases = {
            'b1': tf.Variable(tf.random_normal([num_inputs])),
            'b2': tf.Variable(tf.random_normal([num_inputs])),
            'out': tf.Variable(tf.random_normal([num_classes]))
        }
        # Hidden layer with RELU activation
        layer_1 = tf.add(tf.matmul(x, weights['h1']), biases['b1'])
        layer_1 = tf.nn.sigmoid(layer_1)
        # Hidden layer with RELU activation
        layer_2 = tf.add(tf.matmul(layer_1, weights['h2']), biases['b2'])
        layer_2 = tf.nn.sigmoid(layer_2)
        # Output layer with linear activation
        out_layer = tf.add(tf.matmul(layer_2, weights['out']), biases['out'])
        return out_layer

    def get_previous_job_directory(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        return os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/' + bed_file_name + '/31/MostLikelyBreakPointsJob/'))

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
