import io
import os
import pwd
import sys
import json
import time

from kmer import (
    config,
)

import colorama

print('importing bed.py')
# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

class BedTrack:

    def __init__(self, chrom, begin, end):
        self.chrom = chrom
        self.begin = begin
        self.end = end

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

def track_from_name(name):
    tokens = name.lower().split('_')
    return BedTrack(tokens[0], int(tokens[1]), int(tokens[2]))
