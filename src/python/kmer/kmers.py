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

import socket
import struct

# ============================================================================================================================ #
# FASTQ helpers
# ============================================================================================================================ #

def parse_fastq(path):
    name = None
    #tracemalloc.start()
    HEADER_LINE = 0
    SEQUENCE_LINE = 1
    THIRD_LINE = 2
    QUALITY_LINE = 3
    state = HEADER_LINE
    # need to skip invalid lines
    fastq_file = open(path)
    line = fastq_file.readline().strip()
    ahead = fastq_file.readline().strip()
    n = 0
    m = 0
    t = time.time()
    while ahead:
        if state == HEADER_LINE:
            if line[0] == '@' and ahead[0] != '@':
                name = line[:-1] # ignore the EOL character
                state = SEQUENCE_LINE
        elif state == SEQUENCE_LINE:
            state = THIRD_LINE
            seq = line[:-1] # ignore the EOL character
            yield seq, name
        elif state == THIRD_LINE:
            state = QUALITY_LINE
        elif state == QUALITY_LINE:
            state = HEADER_LINE
        line = ahead
        ahead = fastq_file.readline()

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

def extract_chromosome(chromosomes):
    c = config.Configuration()
    sequence = ''
    ref = open(c.reference_genome)
    count = 0
    line = ref.readline().lower().strip()
    chrom = 'chr' + str(chromosomes.pop(0))
    found = False
    while True:
        count += 1
        if line == '>' + chrom:
            print('found', chrom)
            while True:
                line = ref.readline().lower().strip()
                count += 1
                if line.startswith('>') or len(line) == 0:
                    print(line)
                    yield sequence, chrom
                    sequence = ''
                    if len(chromosomes) == 0:
                        return
                    chrom = 'chr' + str(chromosomes.pop(0))
                    found = True
                    break
                sequence += line
        if found:
            found = False
            continue
        line = ref.readline().lower().strip()

def get_kmer_count(kmer, index, ref):
    start = time.time()
    c = config.Configuration()
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    port = c.count_server_port
    s.connect(('localhost', port + index))
    s.send(bytearray(kmer, 'ascii'))
    response = s.recv(4) # integer size
    count = struct.unpack('!i', response)[0]
    end = time.time()
    print(end - start)
    return count

def get_novel_kmers(kmers, index):
    novel_kmers = {}
    for kmer in kmers:
        if not kmer in novel_kmers:
            count = get_kmer_count(kmer, index, True)
            if count != 0:
                # this is a novel kmer
                novel_kmers[kmer] = True
    return novel_kmers

def has_novel_kmers(kmers, index):
    # checks if this candidate has a kmer that has not occured in the reference genome
    for kmer in kmers:
        if is_kmer_novel(kmer, index):
            return True
    return False

def is_kmer_novel(kmer, index):
    count = get_kmer_count(kmer, index, True)
    return count == 0

def has_unique_novel_kmers(track, candidate, kmers, index):
    # checks if this candidate has a novel kmer that hasn't occurred in any other candidate
    for kmer in kmers:
        if is_kmer_novel(kmer, index):
            found = False
            for break_point in track:
                # skip the wicked candidate count key
                if break_point.find('candidates') != -1:
                    continue
                if candidate != break_point:
                    # this kmer appears in at least one other break point so no need to go further
                    if kmer in track[break_point]['kmers']:
                        found = True
                        break
            # we didn't find this novel kmer anywhere so it should be unique, no need to check others
            if not found:
                return True
    # we haven't returned yet so we didn't find any uniquely novel kmers
    return False

def get_canonical_kmer_representation(kmer):
    kmer = kmer.upper()
    reverse_complement = reverse_complement_sequence(kmer)
    return kmer if kmer < reverse_complement else reverse_complement

def c_extract_canonical_kmers(k, counter = lambda x: 1, count = 1, *args):
    kmers = {}
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = get_canonical_kmer_representation(s[i : i + k])
            c = counter(kmer)
            if c > count:
                continue
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    return kmers

def extract_canonical_kmers(k, *args):
    kmers = {}
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = get_canonical_kmer_representation(s[i : i + k])
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

def c_extract_kmers(ksize, counter = lambda x: 1, count = 1, *args):
    kmers = {}
    for s in args:
        for i in range(0, len(s) - ksize + 1):
            kmer = s[i : i + ksize]
            if counter(kmer) > count:
                continue
            k = find_kmer(kmer, kmers)
            if not k:
                kmers[kmer] = 1
            else:
                kmers[k] += 1
    return kmers

def extract_kmers(ksize, *args):
    kmers = {}
    for s in args:
        for i in range(0, len(s) - ksize + 1):
            kmer = s[i : i + ksize]
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    return kmers

def stream_kmers(k, *args):
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
