
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
import argparse
import operator
import traceback
import deepdiff


import subprocess

from itertools import chain

from nebula import (
    bed,
    config,
    counter,
    reduction,
    map_reduce,
    statistics,
    visualizer,
)

import pysam

from nebula.debug import *
from nebula.kmers import *
from nebula.logger import *
from nebula.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class PacBioVerificationJob(map_reduce.Job):

    _name = 'PacBioVerificationJob'
    _category = 'genotyping'
    _previous_job =  None

    def load_inputs(self):
        c = config.Configuration()
        self.tracks = c.tracks
        self.alignments = pysam.AlignmentFile(c.bam, "rb")
        self.round_robin(self.tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        slack = 120
        t = c.tracks[track_name]
        reads = self.alignments.fetch(track.chrom, track.begin - slack, track.end + slack)
        bases = []
        if track.chrom != 'chr1' or track.svtype != 'DEL' or not 'DEL' in track_name:
            return None
        ref = extract_chromosome(track.chrom)
        for i in range (track.begin, track.end):
            bases.append({'Del': 0, 'Ref': 0, 'Alt': 0})
        for read in reads:
            seq = read.query_sequence
            base = read.reference_start
            ref_index = 0
            seq_index = 0
            for ct in read.cigartuples:
                o = ct[0]
                l = ct[1]
                for i in range(l):
                    if base + ref_index >= track.begin and base + ref_index < track.end:
                        if o == 0 or o == 7 or o == 8: # aligend; match or missmatch
                            a = 'Ref' if seq[seq_index] == ref[base + ref_index] else 'Alt'
                        elif o == 2: # deleted from reference
                            a = 'Del'
                        elif o == 3: # Skipped (N)
                            a = 'Alt'
                        else:
                            a = None
                        j = base + ref_index - track.begin
                        if a in bases[j]:
                            bases[j][a] += 1
                    if o != 6 and o != 5 and o != 4 and o != 1: # no reference consumed
                        ref_index += 1
                    if o != 6 and o != 5 and o != 3 and o != 2: # no query consumed
                        seq_index += 1
        consensus = {'Del': 0, 'Ref': 0, 'Alt': 0}
        for i, base in enumerate(bases):
            for k in consensus:
                consensus[k] += base[k]
            #bases[i] = bases[i].items()
            #bases[i] = sorted(bases[i], key = lambda x: x[1], reverse = True)
            #print(bases[i])
        l = float(len(bases))
        #consensus = {k: v / l for k, v in sorted(consensus.it, key = lambda x: x[1], reverse = True)}
        print(consensus)
        consensus = {k: int(round(math.log10(v))) if v != 0 else -1 for k, v in consensus.items()}
        print(consensus)
        if consensus['Del'] < consensus['Ref'] and consensus['Del'] < consensus['Alt']:
            consensus['genotype'] = '0/0'
        elif consensus['Del'] > consensus['Ref'] and consensus['Del'] > consensus['Alt']:
            consensus['genotype'] = '1/1'
        elif consensus['Del'] == consensus['Ref'] or consensus['Del'] == consensus['Alt']:
            consensus['genotype'] = '1/0'
        else:
            consensus['genotype'] = '1/0'
        if consensus['genotype'].replace('/', '') != t.lp_genotype:
            system_print_error('Wrong')
        return consensus

    def reduce(self):
        c = config.Configuration()
        m = 0
        n = 0
        with open(os.path.join(self.get_current_job_directory(), 'genotyped.bed'), 'w') as bed_file:
            bed_file.write('\t'.join(['#CHROM', 'BEGIN', 'END', 'LP_GENOTYPE', 'PACBIO_GENOTYPE', '\n']))
            for batch in self.load_output():
                for track in batch:
                    t = c.tracks[track]
                    n += 1
                    g = batch[track]['genotype']
                    if t.lp_genotype == g.replace('/', ''):
                        m += 1
                    elif '1' in t.lp_genotype and '1' in g:
                        m += 1
                    bed_file.write('\t'.join([str(s) for s in [t.chrom, t.begin, t.end, t.lp_genotype, g, '\n']]))
        print(m / float(n))




