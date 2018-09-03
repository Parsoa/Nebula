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
import random
import argparse
import traceback
import statistics as stats
import SocketServer as socketserver

from kmer import (
    bed,
    sets,
    config,
)

from kmer.commons import *

import socket
import struct

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

chroms = {}

@Memoize
def extract_chromosome(chromosome):
    if chromosome in chroms:
        return chroms[chromosome]
    c = config.Configuration()
    sequence = ''
    ref = open(c.reference_genome)
    line = ref.readline().lower().strip()
    found = False
    while True:
        if line.startswith('>chr'):
            chrom = line[line.find('>') + 1:]
            if chrom == chromosome:
                print('extracting ' + chrom)
                while True:
                    line = ref.readline().lower().strip()
                    if line.startswith('>') or len(line) == 0:
                        print(line)
                        return sequence
                    sequence += line.upper()
        line = ref.readline().lower().strip()
        if len(line) == 0:
            break

def extract_chromosomes(chromosomes):
    c = config.Configuration()
    sequence = ''
    ref = open(c.reference_genome)
    line = ref.readline().lower().strip()
    found = False
    while True:
        if line.startswith('>chr'):
            chrom = line[line.find('>') + 1:]
            if chrom in chromosomes:
                print('extracting ' + chrom)
                while True:
                    line = ref.readline().lower().strip()
                    if line.startswith('>') or len(line) == 0:
                        print(line)
                        yield sequence, chrom
                        sequence = ''
                        found = True
                        break
                    sequence += line.upper()
        if found:
            found = False
            continue
        line = ref.readline().lower().strip()
        if len(line) == 0:
            break

def extract_whole_genome():
    c = ['chr' + str(x) for x in range(1, 23)]
    c.append('chrx')
    c.append('chry')
    for seq, chrom in extract_chromosomes(c):
        chroms[chrom] = seq
    return chroms

def canonicalize(seq):
    seq = seq.upper()
    reverse_complement = reverse_complement_sequence(seq)
    return seq if seq < reverse_complement else reverse_complement

def c_extract_canonical_kmers(k = 31, counter = lambda x: 1, count = 1, overlap = True, *args):
    c = config.Configuration()
    kmers = {}
    for s in args:
        i = 0
        while i <= len(s) - k and i >= 0:
            kmer = canonicalize(s[i : i + k])
            if counter(kmer) > count:
                i += 1
                continue
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
            if not overlap:
                i += k
            else:
                i += 1
    return kmers

def extract_canonical_kmers(k = 31, *args):
    c = config.Configuration()
    kmers = {}
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = canonicalize(s[i : i + k])
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    return kmers

def extract_canonical_gapped_kmers(*args):
    c = config.Configuration()
    kmers = {}
    for s in args:
        for i in range(0, len(s) - 35 + 1):
            kmer = canonicalize(s[i : i + 35])
            kmer = kmer[:15] + kmer[20:]
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    return kmers

def find_kmer(k, kmers):
    rc = reverse_complement(k)
    if k in kmers:
        return k
    if rc in kmers:
        return rc
    return None

def c_extract_kmers(k = 31, counter = lambda x: 1, count = 1, *args):
    c = config.Configuration()
    kmers = {}
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = s[i : i + k]
            if counter(kmer) > count:
                continue
            k = find_kmer(kmer, kmers)
            if not k:
                kmers[kmer] = 1
            else:
                kmers[k] += 1
    return kmers

def extract_kmers(k = 31, *args):
    c = config.Configuration()
    kmers = {}
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = s[i : i + k]
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    return kmers

def stream_kmers(k = 31, *args):
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = s[i : i + k]
            yield kmer

def reverse_complement(seq):
    return complement_sequence(seq[::-1])

def reverse_complement_sequence(seq):
    return complement_sequence(seq[::-1])

def complement_sequence(seq):
    # A-> C and C->A
    seq = seq.replace('A', 'Z')
    seq = seq.replace('T', 'A')
    seq = seq.replace('Z', 'T')
    #
    seq = seq.replace('G', 'Z')
    seq = seq.replace('C', 'G')
    seq = seq.replace('Z', 'C')
    #
    return seq

def is_subsequence(x, y):
    it = iter(y)
    return all(c in it for c in x)
