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
            self.tracks = self.generate_random_intervals(ref, 1000)
        else:
            self.tracks = self.load_structural_variations()
        self.assign_event_zygosities()
        self.export_bed(self.tracks, 'all')
        self.export_bed(self.absent, 'absent')
        self.export_bed(self.present, 'present')
        self.export_bed(self.homozygous, 'homozygous')
        self.export_bed(self.heterozygous, 'heterozygous')
        #
        self.extract_whole_genome()
        self.round_robin(self.chroms)

    def extract_whole_genome(self):
        self.chroms = extract_whole_genome()
        for chrom, seq in self.chroms.iteritems():
            print(chrom, 'with', len(seq), 'bases')
            with open(os.path.join(self.get_current_job_directory(), chrom + '.fa'), 'w') as chrom_file:
                chrom_file.write('>' + chrom + '\n')
                chrom_file.write(seq)

    def generate_random_intervals(self, n):
        c = config.Configuration()
        random.seed(c.seed)
        print('generating', n, 'random events')
        l = len(self.ref)
        tracks = []
        offset = 100000
        for i in range(0, n):
            begin = random.randint(offset, l - offset)
            end = begin + random.randint(100, 2000)
            track = bed.BedTrack(chrom = c.chrom, begin = begin, end = end)
            tracks.append(track)
        return self.filter_overlapping_intervals(\
                    sorted(sorted(tracks, key = lambda x: x.begin), key = lambda y: y.chrom)\
                )

    def load_structural_variations(self):
        c = config.Configuration()
        print('loading SVs from file', c.bed_file)
        tracks = []
        for track in bed.load_tracks_from_file(c.bed_file, [('zygosity', None)]):
            if track.end - track.begin > 1000000:
                print(red('too large, skipping', track))
                continue
            tracks.append(track)
        print('Total number of tracks:', len(tracks))
        return self.filter_overlapping_intervals(\
                    sorted(sorted(tracks, key = lambda x: x.begin), key = lambda y: y.chrom)\
                )

    def filter_overlapping_intervals(self, intervals):
        remove = []
        i = 0
        while i < len(intervals):
            for j in range(i + 1, len(intervals)):
                # j is contained inside i
                if intervals[j].chrom != intervals[i].chrom:
                    i = j
                    break
                if intervals[j].begin <= intervals[i].end:
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

    def assign_event_zygosities(self):
        c = config.Configuration()
        if c.random:
            n = len(self.tracks) / 2
            r = sorted(random.sample(xrange(len(self.tracks)), n))
            a = list(filter(lambda i: i not in r, range(len(self.tracks))))
            self.absent = [ self.tracks[i] for i in a ]
            self.present = [ self.tracks[i] for i in r ]
            self.homozygous, self.heterozygous = self.select_events(self.present)
        else:
            self.absent = []
            self.present = []
            self.homozygous = []
            self.heterozygous = []
            for track in self.tracks:
                if track.zygosity:
                    if track.zygosity == '1|0' or track.zygosity == '0|1':
                        self.heterozygous.append(track)
                        self.present.append(track)
                    elif track.zygosity == '1|1':
                        self.homozygous.append(track)
                        self.present.append(track)
                    else:
                        self.absent.append(track)
                else:
                    r = random.randint(0, 2)
                    if r == 2:
                        self.absent.append(track)
                    else:
                        self.present.append(track)
                        if r == 1:
                            self.heterozygous.append(track)
                        else:
                            self.homozygous.append(track)
        print(len(self.tracks), 'non-overlapping tracks')
        print(len(self.absent), '0/0')
        print(len(self.present), '1/1 or 1/0')
        print(len(self.homozygous), '1/1')
        print(len(self.heterozygous), '1/0')

    def export_bed(self, intervals, name):
        print('Exporting BED file:', green(name + '.bed'))
        with open(os.path.join(self.get_current_job_directory(), name + '.bed'), 'w') as bed_file:
            for interval in intervals:
                bed_file.write(interval.export()) 

    def transform(self, seq, chrom):
        print('simulating', chrom)
        self.export_chromosome_fasta(chrom, self.present, chrom + '_strand_1')
        self.export_chromosome_fasta(chrom, self.homozygous, chrom + '_strand_2')
        self.export_fastq(seq, chrom + '_strand_1')
        self.export_fastq(seq, chrom + '_strand_2')
        exit()

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

    def apply_events_to_chromosome(self, chrom, tracks):
        seq = ''
        previous = 0
        for track in tracks:
            if track.chrom != chrom:
                print(yellow('skipping', track, 'to', chrom))
                continue
            print('applying', track, 'to', chrom)
            s = self.chroms[chrom][previous:track.begin]
            seq += s
            previous = track.end
        seq += self.chroms[chrom][previous:]
        return seq

    def export_fastq(self, seq, name):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating FASTQ files:', green(name + '.1.fq'), green(name + '.2.fq'))
        num_reads = len(seq) * c.simulation / 100
        fasta = os.path.join(self.get_current_job_directory(), name + '.fa')
        fastq_1 = os.path.join(self.get_current_job_directory(), name + '.1.fq')
        fastq_2 = os.path.join(self.get_current_job_directory(), name + '.2.fq')
        command =  "wgsim -d400 -N{} -1100 -2100 -r{} -R0 -e{} {} {} {}".format(num_reads, c.mutation_rate, c.sequencing_error_rate, fasta, fastq_1, fastq_2)
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def reduce(self):
        self.export_reference_genome()
        self.export_reference_jellyfish_table()
        self.merge_fastq_files()

    def export_reference_genome(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        command = 'cat '
        for chrom in self.chroms:
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
