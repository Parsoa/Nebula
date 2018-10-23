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
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Simulation(map_reduce.Job):

    _name = 'Simulation'
    _category = 'simulation'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = Simulation(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # Simulation
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        if c.seed:
            random.seed(c.seed)
        if c.random:
            self.SVs = self.generate_random_intervals(ref, 1000)
            n = len(self.SVs) / 2
            r = sorted(random.sample(xrange(len(self.SVs)), n))
            a = list(filter(lambda i: i not in r, range(len(self.SVs))))
            self.absent = [ self.SVs[i] for i in a ]
            self.present = [ self.SVs[i] for i in r ]
        else:
            self.SVs = self.load_structural_variations()
        self.export_bed(self.SVs, 'all')
        self.export_bed(self.absent, 'absent')
        self.export_bed(self.present, 'present')
        self.homozygous, self.heterozygous = self.select_events(self.present)
        self.export_bed(self.homozygous, 'homozygous')
        self.export_bed(self.heterozygous, 'heterozygous')
        #
        self.extract_whole_genome()
        self.round_robin(self.chrom)

    def extract_whole_genome(self):
        c = ['chr' + str(x) for x in range(1, 23)]
        c.append('chrx')
        c.append('chry')
        self.chrom = {}
        n = 0
        for seq, chrom in extract_chromosomes(c):
            print('got', chrom, len(seq), 'bases')
            self.chrom[chrom] = seq
            n += 1
            with open(os.path.join(self.get_current_job_directory(), chrom + '.fa'), 'w') as chrom_file:
                chrom_file.write('>' + chrom + '\n')
                chrom_file.write(seq)

    def generate_random_intervals(self, n):
        c = config.Configuration()
        random.seed(c.seed)
        print('generating', n, 'random events')
        l = len(self.ref)
        intervals = []
        offset = 100000
        for i in range(0, n):
            begin = random.randint(offset, l - offset)
            end = begin + random.randint(100, 2000)
            interval = bed.BedTrack(chrom = c.chrom, begin = begin, end = end)
            intervals.append(interval)
        intervals = sorted(intervals, key = lambda x: x.begin)
        intervals = self.filter_overlapping_intervals(intervals)
        return intervals

    def load_structural_variations(self):
        c = config.Configuration()
        print(c.bed_file)
        # split events into batches
        tracks = []
        for track in bed.load_tracks_from_file(c.bed_file, ['zygosity']):
            if track.end - track.begin > 1000000:
                print(red('too large, skipping', track))
                continue
            tracks.append(track)
        tracks = sorted(tracks, key = lambda x: x.begin)
        print('Total number of deletions:', len(tracks))
        return self.filter_overlapping_intervals(tracks)

    def filter_overlapping_intervals(self, intervals):
        remove = []
        i = 0
        while i < len(intervals):
            for j in range(i + 1, len(intervals)):
                # j is contained inside i
                if intervals[j].begin + intervals[j].left_drift <= intervals[i].end + intervals[i].right_drift:
                    remove.append(j)
                    print(red(str(intervals[j])), 'overlaps', blue(str(intervals[i])))
                    continue
                else:
                    i = j
                    break
            if i == len(intervals) - 1:
                break
        n = 0
        for index in sorted(remove):
            intervals.pop(index - n)
            n = n + 1
        return intervals

    def export_bed(self, intervals, name):
        print('Exporting BED file:', green(name + '.bed'))
        with open(os.path.join(self.get_current_job_directory(), name + '.bed'), 'w') as bed_file:
            for interval in intervals:
                bed_file.write(str(interval.chrom) + '\t' + str(interval.begin) + '\t' + str(interval.end) + '\t' +
                    str(interval.left_drift) + '\t' + str(interval.right_drift) + '\n')

    def select_events(self, SVs):
        c = config.Configuration()
        print('Present SVs:', blue(len(SVs)))
        # select a set of commons SVs to be applied to both strands
        n = len(SVs) / 2
        print('Homozygous SVs:', blue(n))
        SVs = sorted(SVs, key = lambda x: x.begin)
        homozygous = [ SVs[i] for i in sorted(random.sample(range(len(SVs)), n)) ]
        heterozygous = []
        for sv in SVs:
            if not sv in homozygous:
                heterozygous.append(sv)
        # sort
        homozygous = sorted(homozygous, key = lambda x: x.begin)
        heterozygous = sorted(heterozygous, key = lambda x: x.begin)
        return homozygous, heterozygous

    def transform(self, seq, chrom):
        print('simulating', chrom)
        self.export_chromosome_fasta(chrom, self.present, chrom + '_strand_1')
        self.export_chromosome_fasta(chrom, self.homozygous, chrom + '_strand_2')
        self.export_fastq(seq, chrom + '_strand_1')
        self.export_fastq(seq, chrom + '_strand_2')

    def export_chromosome_fasta(self, chrom, events, name):
        print('Exporting FASTA file:', green(name + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), name + '.fa'), 'w') as fasta_file:
            fasta_file.write('>' + chrom + '\n')
            seq = self.apply_events_to_chromosome(chrom, events)
            l = len(seq)
            n = 100
            num_lines = l / n
            for i in range(0, num_lines):
                line = seq[i * n : (i + 1) * n].upper() + '\n'
                fasta_file.write(line)

    def apply_events_to_chromosome(self, chrom, SVs):
        seq = ''
        real = 0
        previous = 0
        index = {}
        for i, sv in enumerate(SVs):
            left = sv.begin + sv.left_drift
            s = self.chrom[chrom][previous:left]
            seq += s
            previous = sv.end + sv.right_drift
        seq += chrom[previous:]
        return seq

    def export_fastq(self, seq, name):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating FASTQ files:', green(name + '.1.fq'), green(name + '.2.fq'))
        num_reads = len(seq) * c.simulation / 100
        fasta = os.path.join(self.get_current_job_directory(), name + '.fa')
        fastq_1 = os.path.join(self.get_current_job_directory(), name + '.1.fq')
        fastq_2 = os.path.join(self.get_current_job_directory(), name + '.2.fq')
        command =  "wgsim -d400 -N{} -1100 -2100 -r0 -R0 -e0.001 {} {} {}".format(num_reads, fasta, fastq_1, fastq_2)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def reduce(self):
        self.chrom = ['chr' + c for c in [str(x) for x in range(1, 23)]]
        self.chrom.append('chrx')
        self.chrom.append('chry')
        self.export_reference_genome()
        self.export_reference_jellyfish_table()
        self.merge_fastq_files()

    def export_reference_genome(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        command = 'cat '# + os.path.join(self.get_current_job_directory(), 'chr*.fq') + ' > ' + os.path.join(self.get_current_job_directory(), 'test.fq')
        for chrom in self.chrom:
            path = os.path.join(self.get_current_job_directory(), chrom + '.fa')
            command += path + ' ' 
        path = os.path.join(self.get_current_job_directory(), 'reference.fa')
        command += '> ' + path
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def export_reference_jellyfish_table(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating Jellyfish table')
        command = "jellyfish count -m " + str(c.ksize) + " -s 1000000000 -t 24 --canonical --out-counter-len 2 "
        command += os.path.join(self.get_current_job_directory(), 'reference.fa')
        command += ' -o ' + os.path.join(self.get_current_job_directory(), 'reference_' + str(c.ksize) + '.jf')
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def merge_fastq_files(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Mixing FASTQ files...')
        command = 'cat '
        command += os.path.join(self.get_current_job_directory(), 'chr*.1.fq')
        command += ' > '
        command += os.path.join(self.get_current_job_directory(), 'test.fq')
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def export_fastq_jellyfish_table(self, channel, *args):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating Jellyfish table')
        command = "jellyfish count -m 31 -s 1000000000 -t 24 --canonical --out-counter-len 2"
        for name in args:
            command += ' ' + os.path.join(self.get_current_job_directory(), name + '.fq')
        command += ' -o ' + os.path.join(self.get_current_job_directory(), channel + '.jf')
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
