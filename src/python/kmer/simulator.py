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

from kmer import (
    bed,
    sets,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

def import_reference_genome():
    c = config.Configuration()
    sequence = ''
    ref = open(c.reference_genome)
    count = 0
    print("Copying chromesome {} from {}".format(c.chrom, c.reference_genome))
    while True:
        count += 1
        line = ref.readline().strip()
        if line == '>' + c.chrom or line == '':
            break
    if len(line) == 0:
        return []
    print('Found chromosome', c.chrom, 'at line', count)
    while True:
        count += 1
        if count % 10000000 == 0:
            print("\tReading line {} ...".format(count))
        line = ref.readline().upper().strip()
        if line.startswith('>') or line == '':
            break
        sequence += line
    print("Completed. There are {} bases in chromosome {}.".format(len(sequence), c.chrom))
    return sequence

def load_structural_variations():
    c = config.Configuration()
    bedtools = pybedtools.BedTool(c.bed_file)
    # split events into batches
    n = 0
    tracks = []
    for track in bedtools:
        if track.chrom != c.chrom:
            continue
        name = re.sub(r'\s+', '_', str(track).strip()).strip()
        # too large, skip
        if track.end - track.start > 1000000:
            print(red('skipping', name))
            continue
        tracks.append(track)
    tracks = sorted(tracks, key = lambda x: x.start)
    print('Total number of deletions:', len(tracks))
    return filter_overlapping_intervals(tracks)

def filter_overlapping_intervals(intervals):
    remove = []
    i = 0
    while i < len(intervals):
        for j in range(i + 1, len(intervals)):
            # j is contained inside i
            if intervals[j].end < intervals[i].end:
                remove.append(j)
                print(red(str(intervals[j])), 'overlaps', blue(str(intervals[i])))
                continue
            if intervals[j].start > intervals[i].end:
                i = j
                break
        if i == len(intervals) - 1:
            break
    n = 0
    for index in sorted(remove):
        intervals.pop(index - n)
        n = n + 1
    return intervals

def generate_fasta(ref, SVs):
    c = config.Configuration()
    print('Total number of SVs:', blue(len(SVs)))
    # sample from these SVs
    #n = len(SVs) / 5
    #print('Sampling', blue(n), 'events')
    #SVs = [ SVs[i] for i in sorted(random.sample(xrange(len(SVs)), n)) ]
    # select a set of commons SVs to be applied to both
    n = len(SVs) / 2
    print('Shared events:', blue(n))
    common = [ SVs[i] for i in sorted(random.sample(xrange(len(SVs)), n)) ]
    # sort (just in case)
    SVs = sorted(SVs, key = lambda x: x.start)
    common = sorted(common, key = lambda x: x.start)
    test = ''
    control = ''
    ref_size = len(ref)
    # apply all events to control sequence
    previous = 0
    for i, sv in enumerate(SVs):
        right = sv.end
        left = sv.start
        control += ref[previous:left]
        previous = right
    # only apply commons SVs to test sequence
    previous = 0
    for i, sv in enumerate(common):
        right = sv.end
        left = sv.start
        test += ref[previous:left]
        previous = right
    return control, test, SVs, common

def export_fasta(sequence, name):
    print('Exporting FASTQ file:', green(name)) 
    c = config.Configuration()
    for i in range(0, 2):
        with open(os.path.join(get_output_directory(), name + '.fa'), 'w') as fasta_file:
            fasta_file.write('>' + c.chrom + '_' + str(i) + '\n')
            l = len(sequence)
            num_lines = l / 50
            for i in range(0, num_lines):
                fasta_file.write(sequence[i * 50 : (i + 1) * 50] + '\n')
            if l % 50 != 0:
                fasta_file.write(sequence[num_lines * 50 :] + '\n')
    print('Done!')

def export_bam(name):
    print('Exporting BAM file:', green(name + '.bam'))
    c = config.Configuration()
    # gerenate fq files
    print('Generating FASTQ file ...')
    num_reads = 1000000 * c.coverage
    fasta = os.path.join(get_output_directory(), name + '.fa')
    fastq_1 = os.path.join(get_output_directory(), name + '.1.fq')
    fastq_2 = os.path.join(get_output_directory(), name + '.2.fq')
    command =  "wgsim -d400 -N{} -1100 -2100 {} {} {}".format(num_reads, fasta, fastq_1, fastq_2)
    os.system(command)
    # generate sam file
    print('Generating SAM file ...')
    sam_file = os.path.join(get_output_directory(), name + '.sam')
    command = "bwa mem -M -t 20 -R \"@RG\\tID:1\\tPL:ILLUMINA\\tSM:cnv_1000_ref\" {} {} {} > {}".format(c.reference_genome, fastq_1, fastq_2, sam_file)
    os.system(command)
    # generate bam file
    prina('Generating unsorted BAM file ...')
    unsorted_bam = os.path.join(get_output_directory(), name + '.unsorted.bam')
    command = "samtools view -S -b {} > {}".format(sam_file, unsorted_bam)  
    os.system(command)
    print('Sorting ...')
    bam = os.path.join(get_output_directory(), name + '.bam')
    command = "samtools sort {} {}".format(unsorted_bam, bam)
    os.system(command)
    print('Indexing ...')
    bam_index = os.path.join(get_output_directory(), name + '.bai')
    command = "samtools index {}".format(bam_index)
    os.system(command)
    print('Done!')

def get_output_directory():
    return os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/simulation/'))

def export_bed(intervals, name):
    print('Exporting BED file:', green(name + '.bed'))
    with open(os.path.join(get_output_directory(), name + '.bed'), 'w') as bed_file:
        for interval in intervals:
            bed_file.write(str(interval.chrom) + '\t' + str(interval.start) + '\t' + str(interval.end) + '\n')

def simulate():
    c = config.Configuration()
    c.chrom = 'chr1'
    ref = import_reference_genome()
    SVs = load_structural_variations()
    control, test, SVs, present = generate_fasta(ref, SVs)
    export_bed(SVs, 'all')
    export_bed(present, 'present')
    export_fasta(control, 'control')
    export_fasta(test, 'test')
    export_bam('control')
    export_bam('test')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    start = time.time()
    simulate()
    end = time.time()
    print("Simulation time: {} seconds".format(end-start))
