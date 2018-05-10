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
import operator
import traceback
import subprocess

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

def extract_chromosome(chrom):
    c = config.Configuration()
    sequence = ''
    ref = open(c.reference_genome)
    count = 0
    while True:
        count += 1
        line = ref.readline().strip()
        if line == '>' + c.chrom or line == '':
            break
    if len(line) == 0:
        return []
    while True:
        count += 1
        line = ref.readline().upper().strip()
        if line.startswith('>') or line == '':
            break
        sequence += line
    return sequence

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

def generate_random_intervals(ref, n):
    c = config.Configuration()
    random.seed(c.seed)
    print('generating', n, 'random events')
    l = len(ref)
    SVs = []
    offset = 100000
    for i in range(0, n):
        start = random.randint(offset, l - offset)
        end = start + random.randint(100, 2000)
        event = pybedtools.Interval(chrom = 'chr1', start = start, end = end)
        SVs.append(event)
    SVs = sorted(SVs, key = lambda x: x.start)
    SVs = filter_overlapping_intervals(SVs)
    return SVs

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
            if intervals[j].start >= intervals[i].start and intervals[j].start <= intervals[i].end:
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
    # select a set of commons SVs to be applied to both
    n = len(SVs) / 2
    print('Present SVs:', blue(n))
    common = [ SVs[i] for i in sorted(random.sample(xrange(len(SVs)), n)) ]
    # sort (just in case)
    SVs = sorted(SVs, key = lambda x: x.start)
    common = sorted(common, key = lambda x: x.start)
    test = apply_events_to_ref(ref, common)
    control = apply_events_to_ref(ref, SVs)
    return control, test, SVs, common

def generate_homozygous_fasta(ref, SVs):
    c = config.Configuration()
    print('Total number of SVs:', blue(len(SVs)))
    # sort (just in case)
    SVs = sorted(SVs, key = lambda x: x.start)
    sequence = apply_events_to_ref(ref, SVs)
    return sequence

def generate_heterozygous_fasta(ref, SVs):
    c = config.Configuration()
    print('Present SVs:', blue(len(SVs)))
    # select a set of commons SVs to be applied to both strands
    n = len(SVs) / 2
    print('Homozygous SVs:', blue(n))
    homozygous = [ SVs[i] for i in sorted(random.sample(xrange(len(SVs)), n)) ]
    heterozygous = []
    for sv in SVs:
        if not sv in homozygous:
            heterozygous.append(sv)
    # sort
    SVs = sorted(SVs, key = lambda x: x.start)
    homozygous = sorted(homozygous, key = lambda x: x.start)
    heterozygous = sorted(heterozygous, key = lambda x: x.start)
    strand_1 = apply_events_to_ref(ref, SVs)
    strand_2 = apply_events_to_ref(ref, homozygous)
    return strand_1, strand_2, homozygous, heterozygous

def export_fasta(sequence, name):
    print('Exporting FASTA file:', green(name + '.fa')) 
    c = config.Configuration()
    with open(os.path.join(get_output_directory(), name + '.fa'), 'w') as fasta_file:
        fasta_file.write('>' + c.chrom + '\n')
        l = len(sequence)
        num_lines = l / 50
        for i in range(0, num_lines):
            fasta_file.write(sequence[i * 50 : (i + 1) * 50] + '\n')
        if l % 50 != 0:
            fasta_file.write(sequence[num_lines * 50 :] + '\n')

def get_sv_type():
    c = config.Configuration()
    bed_file_name = c.bed_file.split('/')[-1]
    sv_type = bed_file_name.split('.')[-2]
    return sv_type

def apply_events_to_ref(ref, SVs):
    seq = ''
    previous = 0
    for i, sv in enumerate(SVs):
        right = sv.end
        left = sv.start
        seq += ref[previous:left]
        if get_sv_type() == 'INV':
            seq += reverse_complement_sequence(ref[left:right])
        previous = right
    return seq

def export_bam(name):
    FNULL = open(os.devnull, 'w')
    c = config.Configuration()
    print('Generating FASTQ files:', green(name + '.1.fq'), green(name + '.2.fq'))
    num_reads = 1000000 * c.coverage
    fasta = os.path.join(get_output_directory(), name + '.fa')
    fastq_1 = os.path.join(get_output_directory(), name + '.1.fq')
    fastq_2 = os.path.join(get_output_directory(), name + '.2.fq')
    command =  "wgsim -d400 -N{} -1100 -2100 {} {} {}".format(num_reads, fasta, fastq_1, fastq_2)
    output = subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
    # generate sam file
    print('Generating SAM file:', green(name + '.sam'))
    sam_file = os.path.join(get_output_directory(), name + '.sam')
    command = "bwa mem -M -t 20 -R \"@RG\\tID:1\\tPL:ILLUMINA\\tSM:cnv_1000_ref\" {} {} {} > {}".format(c.reference_genome, fastq_1, fastq_2, sam_file)
    subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
    # generate bam file
    print('Generating unsorted BAM file:', green(name + '.unsorted.bam'))
    unsorted_bam = os.path.join(get_output_directory(), name + '.unsorted.bam')
    command = "samtools view -S -b {} > {}".format(sam_file, unsorted_bam)  
    subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
    print('Sorting ...')
    bam = os.path.join(get_output_directory(), name)
    command = "samtools sort {} -o {}".format(unsorted_bam, bam)
    subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
    print('Indexing ...')
    bam_index = os.path.join(get_output_directory(), name + '.bai')
    command = "samtools index {} {}".format(bam + '.bam', bam_index)
    subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
    print('Done!')

def get_output_directory():
    c = config.Configuration()
    bed_file_name = c.bed_file.split('/')[-1]
    return os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/' + bed_file_name + '/' + str(c.ksize) + '/Simulation'))

def create_output_directories():
    dir = get_output_directory()
    if not os.path.exists(dir):
        os.makedirs(dir)

def export_bed(intervals, name):
    print('Exporting BED file:', green(name + '.bed'))
    with open(os.path.join(get_output_directory(), name + '.bed'), 'w') as bed_file:
        for interval in intervals:
            bed_file.write(str(interval.chrom) + '\t' + str(interval.start) + '\t' + str(interval.end) + '\n')

def simulate():
    c = config.Configuration()
    c.chrom = 'chr1'
    ref = import_reference_genome()
    if c.seed:
        random.seed(c.seed)
    if c.random:
        SVs = generate_random_intervals(ref, 1000)
    else:
        SVs = load_structural_variations()
    if c.heterozygous:
        pid = os.fork()
        # export homozygous control genome
        if pid == 0:
            export_bed(SVs, 'all')
            control = generate_homozygous_fasta(ref, SVs)
            export_fasta(control, 'control')
            export_bam('control')
        # use a second thread to export the test sample
        else:
            n = len(SVs) / 2
            present = [ SVs[i] for i in sorted(random.sample(xrange(len(SVs)), n)) ]
            export_bed(present, 'present')
            strand_1, strand_2, homozygous_SVs, heterozygous_SVs = generate_heterozygous_fasta(ref, present)
            export_bed(homozygous_SVs, 'homozygous')
            export_bed(heterozygous_SVs, 'heterozygous')
            export_fasta(strand_1, 'test_strand_1')
            export_fasta(strand_2, 'test_strand_2')
            # export a couple of fastq files for each strand, will this work?
            export_bam('test_strand_1')
            export_bam('test_strand_2')
    else:
        control, test, SVs, present = generate_fasta(ref, SVs)
        pid = os.fork()
        if pid == 0:
            export_bed(present, 'present')
            export_fasta(test, 'test')
            export_bam('test')
        else:
            export_bed(SVs, 'all')
            export_fasta(control, 'control')
            export_bam('control')

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    start = time.time()
    simulate()
    end = time.time()
    print("Simulation time: {} seconds".format(end-start))
