from __future__ import print_function

import io
import os
import pwd
import sys
import json
import time
import ctypes

from kmer import (
    config,
)

from kmer.kmers import *
from kmer.commons import *

#import khmer
import jellyfish

# ============================================================================================================================ #
# Caching
# ============================================================================================================================ #

#def cache(f):
#    cache = pylru.lrucache(config.Configuration.kmer_cache_size)
#    hits = 0
#    misses = 0
#    def wrapper(kmer, index):
#        if kmer in cache:
#            print('miss')
#            nonlocal misses
#            misses += 1
#            return cache[kmer]
#        nonlocal hits
#        hits += 1
#        print('hit: ', hits / (hits + misses), 'hits: ', hits, ' misses: ', misses)
#        cache[kmer] = f(kmer, index) 
#        return cache[kmer]
#    return wrapper

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class KmerCountsProvider(object):

    def __init__(self):
        self.import_counts()

# ============================================================================================================================ #
# Dummy. Always returns 40. Use for testing.
# ============================================================================================================================ #

class DictionaryCountsProvider(KmerCountsProvider):

    def __init__(self, kmers):
        self.kmers = kmers

    def get_kmer_count(self, kmer):
        k = find_kmer(kmer, self.kmers)
        return self.kmers[k]['count'] if k else None

    def get_kmer(self, kmer):
        k = find_kmer(kmer, self.kmers)
        if k:
            return self.kmers[k]
        return None

    def has_kmer(self, kmer):
        return find_kmer(kmer, self.kmers)

    def stream_kmers(self):
        for kmer in self.kmers:
            yield str(kmer), self.kmers[kmer]['count']

# ============================================================================================================================ #
# Dummy. Always returns 40. Use for testing.
# ============================================================================================================================ #

class ReferenceCountsProvider(KmerCountsProvider):

    def __init__(self, path):
        self.path = path
        self.chrom = {}
        for c in range(1, 23):
            print(c)
            self.chrom[c] = open(os.path.join(self.path, 'chr' + str(c) + '.fa')).readlines()[1].upper()
        self.chrom['x'] = open(os.path.join(self.path, 'chrx.fa')).readlines()[1].upper()
        self.chrom['y'] = open(os.path.join(self.path, 'chry.fa')).readlines()[1].upper()

    def get_kmer_count(self, kmer):
        return sum(list(map(lambda x: len(x), [[m.start() for m in re.finditer(kmer, self.chrom[c])] for c in self.chrom]))) + sum(list(map(lambda x: len(x), [[m.start() for m in re.finditer(reverse_complement(kmer), self.chrom[c])] for c in self.chrom])))

# ============================================================================================================================ #
# Dummy. Always returns 40. Use for testing.
# ============================================================================================================================ #

class DummyCountsProvider(KmerCountsProvider):

    def get_kmer_count(self, kmer):
        return 40

# ============================================================================================================================ #
# Jellyfish
# ============================================================================================================================ #

class JellyfishCountsProvider(KmerCountsProvider):

    def __init__(self, path):
        self.path = path
        self.import_counts()

    def import_counts(self):
        print('importing jellyfish table', self.path)
        self.qf = jellyfish.QueryMerFile(self.path)
        print('table loaded')

    def get_kmer_count(self, kmer):
        canon = jellyfish.MerDNA(str(kmer))
        canon.canonicalize()
        return self.qf[canon]

    def stream_kmers(self):
        mf = jellyfish.ReadMerFile(self.path)
        #print({key: value for key, value in mf.__dict__.items() if not key.startswith("__")})
        for kmer, count in mf:
            yield str(kmer), count

# ============================================================================================================================ #
# KMC
# ============================================================================================================================ #

class KmcCountsProvider(KmerCountsProvider):

    def __init__(self, path):
        self.path = path
        self.import_counts()

    def import_counts(self):
        print('importing jellyfish table')
        self.kmc = ctypes.CDLL(os.path.join(os.path.dirname(__file__), "../../cpp/kmc.so"))
        print('table loaded')

    def get_kmer_count(self, kmer):
        kmer_str = ctypes.create_string_buffer(str.encode(kmer))
        count = self.kmc.get_kmer_count(kmer_str)
        return count

