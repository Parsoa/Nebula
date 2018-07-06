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

class DummyCountsProvider(KmerCountsProvider):

    def get_kmer_count(self, kmer):
        return 40

# ============================================================================================================================ #
# Khmer
# ============================================================================================================================ #

#class KhmerCountsProvider(KmerCountsProvider):
#
#    @commons.measure_time
#    def export_counts():
#        c = config.Configuration()
#        # 
#        cache = c.counttable + '.ct'
#        print(colorama.Fore.BLUE + 'searching for cached counttable ', cache)
#        if os.path.isfile(cache):
#            print(colorama.Fore.BLUE + 'found at ', cache)
#            return
#        #
#        print(colorama.Fore.BLUE + 'not found, generating counttable...')
#        counttable, nkmers = self.count_kmers_from_file(c.khmer)
#        counttable.save(cache)
#        print(colorama.Fore.BLUE + 'done')
#        #
#        print(colorama.Fore.BLUE + 'counttable cached\n', 'kmers: ', nkmers,\
#            '\nsize: ', os.stat(cache).st_size)
#        return
#
#    def count_kmers_from_file(seq_file):
#        c = config.Configuration()
#        #
#        counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
#        nseqs, nkmers = counttable.consume_seqfile(seq_file)
#        #
#        return counttable, nkmers
#
#    @commons.measure_time
#    def import_counts(self):
#        c = config.Configuration()
#        print(colorama.Fore.MAGENTA + 'importing counttable for ', c.counttable)
#        cache = c.counttable + '.ct'
#        counttable = khmer.Counttable.load(cache)
#        print(colorama.Fore.MAGENTA + 'done')
#        return counttable

# ============================================================================================================================ #
# Jellyfish
# ============================================================================================================================ #

class JellyfishCountsProvider(KmerCountsProvider):

    def __init__(self, path):
        self.path = path
        self.import_counts()

    def import_counts(self):
        print('importing jellyfish table')
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
            yield kmer, count

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

