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

chroms = {}

def extract_chromosome(chromosome):
    chromosome = chromosome.lower()
    if chromosome in chroms:
        print(yellow('loading from cache'))
        return chroms[chromosome]
    else:
        print(red('chromosome not found'), chromosome)
    c = config.Configuration()
    sequence = ''
    print(yellow(c.reference_genome))
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
                        chroms[chromosome] = sequence
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
    m = 0
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
                        m += 1
                        if m == len(chromosomes):
                            return
                        break
                    sequence += line.upper()
        # this is to avoid skipping the last line we read for the previous chromosome (header of next)
        if found:
            found = False
            continue
        line = ref.readline().lower().strip()
        if len(line) == 0:
            break

def extract_whole_genome():
    print('extracting whole genome')
    c = ['chr' + str(x) for x in range(1, 23)]
    c.append('chrx')
    c.append('chry')
    for seq, chrom in extract_chromosomes(c):
        chroms[chrom] = seq
    return chroms

