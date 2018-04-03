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
# kmer helpers
# ============================================================================================================================ #

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

def extract_kmers(k, *args):
    kmers = {}
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = s[i : i + k]
            canon = kmer#get_canonical_kmer_representation(kmer)
            if not canon in kmers:
                kmers[canon] = 0
            kmers[canon] += 1
    return kmers

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
