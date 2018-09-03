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
    counter,
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import acora
import cplex
import numpy
import pybedtools

from Bio import pairwise2

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class FindReadsForKmerJob(counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = FindReadsForKmerJob(job_name = 'FindReadsForKmerJob_', previous_job_name = '', category = 'debug', **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.kmers = [
            ''
        ]

    # So either a locus has unique indicator kmers or it doesn't
    # if it doesn't then it means some or all of it's kmers are shared with other loci of the same kmer
    # A target either has unique kmers so the first loop should cover it or it doesn't so the second loop should
    # Read errors and SNPs when a target has unique kmers will also be handled by the second loop
    def transform(self):
        c = config.Configuration()
        #self.fastq_file = open(track, 'r')
        for read, name in self.parse_fastq():
            kmers = extract_canonical_kmers(c.ksize, read)
            found = False
            for kmer in kmers:
                if kmer in self.kmers:
                    self.kmers[kmer]['total'] += 1
                    for ukmer in kmers:
                        if ukmer in self.kmers[kmer]['indicator_kmers']:
                            self.kmers[kmer]['count'] += 1
                            found = True
                            break
                    i = read.find(kmer)
                    if i == -1:
                        i = read.find(reverse_complement(kmer))
                    left = read[:i]
                    right = read[i + c.ksize:]
                    for locus in self.kmers[kmer]['loci']:
                        if found:
                            break
                        for mask in self.kmers[kmer]['loci'][locus]['mask']:
                            if is_subsequence(mask, left) or is_subsequence(mask, right):
                                self.kmers[kmer]['doubt'] += 1
                                found = True
                                break
                                        #self.kmers[kmer]['loci'][locus]['count'] += 1
                    #if self.kmers[kmer]['unique']:
                    #    self.kmers[kmer]['count'] += 1
                    #    continue
                    #n = list(filter(lambda x: x in kmers, self.kmers[kmer]['negative']))
                    #p = list(filter(lambda x: x in kmers, self.kmers[kmer]['positive']))
                    #if p and not n:
                    #    #print(green('positive'))
                    #    self.kmers[kmer]['count'] += 1
                    #    #self.kmers[kmer]['reads']['Pn'].append(read)
                    #if n and not p:
                    #    #print(red('negative'))
                    #    #self.kmers[kmer]['reads']['pN'].append(read)
                    #    pass
                    #else:
                    #    #if p and n:
                    #    #    self.kmers[kmer]['reads']['PN'].append(read)
                    #    #else:
                    #    #    self.kmers[kmer]['reads']['pn'].append(read)
                    #    self.kmers[kmer]['doubt'] += 1
                    #    continue
                    #    l = {}
                    #    for o in self.kmers[kmer]['ocurrences']:
                    #        l[o] = len(list(filter(lambda x: x in kmers, self.kmers[kmer]['ocurrences'][o]['kmers'])))
                    #    if sum(map(lambda o: l[o], l)):
                    #        m = 0
                    #        t = bed.track_from_name(self.kmers[kmer]['track'])
                    #        for position in l:
                    #            tokens = position.split('_')
                    #            if tokens[0] == t.chrom and int(tokens[1]) >= t.start and int(tokens[1]) < t.end:
                    #                m += l[position]
                    #        self.kmers[kmer]['count'] += float(m) / sum(map(lambda o: l[o], l))
                    #    else:
                    #        self.kmers[kmer]['count'] += 1.0 / float(len(l)) 
                    #else:
                    #    self.kmers[kmer]['doubt'] += 1
                    #    #print(yellow('both or none'))
                    #    seq = read[: i] + read[i + c.ksize:]
                    #    choice = (-10000000, [])
                    #    for position in self.kmers[kmer]['ocurrences']:
                    #        ref = self.kmers[kmer]['ocurrences'][position]['seq']['left'] + self.kmers[kmer]['ocurrences'][position]['seq']['right']
                    #        #print(seq)
                    #        #print(reverse_complement(seq))
                    #        #print(ref)
                    #        alignments = pairwise2.align.globalxs(seq, ref, -1, -1)
                    #        score = alignments[0][2]
                    #        choice = (score, [(position, alignments[0])]) if score > choice[0] else (score, choice[1] + [(position, alignments[0])]) if score == choice[0] else choice
                    #        alignments = pairwise2.align.globalxs(reverse_complement(seq), ref, -1, -1)
                    #        score = alignments[0][2]
                    #        choice = (score, [(position, alignments[0])]) if score > choice[0] else (score, choice[1] + [(position, alignments[0])]) if score == choice[0] else choice
                    #    if choice[0] > len(seq) - 10:
                    #        m = 0
                    #        t = self.kmers[kmer]['track']
                    #        for (position, a) in choice[1]:
                    #            tokens = position.split('_')
                    #            if tokens[0] == t.chrom and int(tokens[1]) >= t.start and int(tokens[1]) < t.end:
                    #                m += 1
                    #                #print(pairwise2.format_alignment(*a))
                    #        self.kmers[kmer]['count'] += 1.0 * (float(m) / len(choice[1]))
                    #    # read is too different from any of the kmer's occurrences, assume error
                    #    else:
                    #        self.kmers[kmer]['count'] += 1.0 / self.kmers[kmer]['reference']

    def reduce(self):
        if not self.resume_from_reduce:
            self.kmers = self.merge_counts('count', 'total')
        else:
            self.kmers = json.load(open(os.path.join(self.get_current_job_directory(), 'kmers.json')))
        x = []
        y = []
        for kmer in self.kmers:
            t = len(self.kmers[kmer]['loci'])
            x.append(self.kmers[kmer]['count'] / t)
            y.append(t)
        visualizer.histogram(x = x, name = 'average coverage after reduction', path = self.get_current_job_directory(), x_label = 'average coverage', y_label = 'number of kmers')
        visualizer.histogram(x = y, name = 'loci of interest', path = self.get_current_job_directory(), x_label = 'number of loci', y_label = 'number of kmers')
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

