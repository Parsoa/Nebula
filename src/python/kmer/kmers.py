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

from kmer import (
    config,
)

from kmer.commons import *

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

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

def stream_canonical_kmers(k = 31, *args):
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = s[i : i + k]
            yield canonicalize(kmer)

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
