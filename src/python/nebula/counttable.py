from __future__ import print_function

from nebula import (
    config,
)

from nebula.kmers import *
from nebula.logger import *

import jellyfish

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class KmerCountsProvider(object):

    def __init__(self):
        self.import_counts()

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
        for kmer, count in mf:
            yield str(kmer), count
