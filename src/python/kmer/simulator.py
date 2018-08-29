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
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from kmer.kmers import *
from kmer.commons import *
print = pretty_print

import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Simulation(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = Simulation(job_name = 'Simulation_', previous_job_name = '', category = 'simulation', **kwargs)
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
        else:
            self.SVs = self.load_structural_variations()
        self.export_bed(self.SVs, 'all')
        n = len(self.SVs) / 2
        r = sorted(random.sample(xrange(len(self.SVs)), n))
        a = list(filter(lambda i: i not in r, range(len(self.SVs))))
        self.absent = [ self.SVs[i] for i in a ]
        self.present = [ self.SVs[i] for i in r ]
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
            start = random.randint(offset, l - offset)
            end = start + random.randint(100, 2000)
            interval = pybedtools.Interval(chrom = c.chrom, start = start, end = end)
            intervals.append(interval)
        intervals = sorted(intervals, key = lambda x: x.start)
        intervals = self.filter_overlapping_intervals(intervals)
        return intervals

    def load_structural_variations(self):
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
        return self.filter_overlapping_intervals(tracks)

    def filter_overlapping_intervals(self, intervals):
        remove = []
        i = 0
        while i < len(intervals):
            for j in range(i + 1, len(intervals)):
                # j is contained inside i
                if intervals[j].start <= intervals[i].end:
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
                bed_file.write(str(interval.chrom) + '\t' + str(interval.start) + '\t' + str(interval.end) + '\n')

    def select_events(self, SVs):
        c = config.Configuration()
        print('Present SVs:', blue(len(SVs)))
        # select a set of commons SVs to be applied to both strands
        n = len(SVs) / 2
        print('Homozygous SVs:', blue(n))
        SVs = sorted(SVs, key = lambda x: x.start)
        homozygous = [ SVs[i] for i in sorted(random.sample(range(len(SVs)), n)) ]
        heterozygous = []
        for sv in SVs:
            if not sv in homozygous:
                heterozygous.append(sv)
        # sort
        homozygous = sorted(homozygous, key = lambda x: x.start)
        heterozygous = sorted(heterozygous, key = lambda x: x.start)
        return homozygous, heterozygous

    def transform(self, seq, chrom):
        print('simulating', chrom)
        self.export_chromosome_fasta(chrom, self.present, chrom + '_strand_1')
        self.export_chromosome_fasta(chrom, self.homozygous, chrom + '_strand_2')
        # export a couple of fastq files for each strand, will this work?
        self.export_fastq(seq, chrom + '_strand_1')
        self.export_fastq(seq, chrom + '_strand_2')
        #self.merge_fastq_files(chrom, chrom + '_strand_1.1', chrom + '_strand_2.2')
        #self.merge_fastq_files(chrom + '_strand_2')

    def export_chromosome_fasta(self, chrom, events, name):
        print('Exporting FASTA file:', green(name + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), name + '.fa'), 'w') as fasta_file:
            fasta_file.write('>' + chrom + '\n')
            if chrom == c.chrom:
                seq = self.apply_events_to_chromosome(chrom, events)
            else:
                seq = self.chrom[chrom]
            l = len(seq)
            n = 100
            num_lines = l / n
            for i in range(0, num_lines):
                line = seq[i * n : (i + 1) * n].upper() + '\n'
                fasta_file.write(line)

    def apply_events_to_chromosome(self, chrom, SVs):
        seq = ''
        previous = 0
        for i, sv in enumerate(SVs):
            left = sv.start
            #print(previous, '\t', left)
            s = self.chrom[chrom][previous:left]
            seq += s
            previous = sv.end
        #print(previous)
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
        command =  "wgsim -d400 -N{} -1100 -2100 -r0 -R0 -e0 {} {} {}".format(num_reads, fasta, fastq_1, fastq_2)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def merge_fastq_files(self, name, *args):
        print('Merging fastq files', name)
        with open(os.path.join(self.get_current_job_directory(), name + '.fq'), 'w') as out_file:
            for arg in args:
                with open(os.path.join(self.get_current_job_directory(), arg + '.fq'), 'r') as in_file:
                    line = in_file.readline()
                    while line:
                        out_file.write(line)
                        line = in_file.readline()
        #os.rename(os.path.join(self.get_current_job_directory(), name + '.1.fq'), os.path.join(self.get_current_job_directory(), name + '.fq'))
        #os.remove(os.path.join(self.get_current_job_directory(), name + '.2.fq'))

    def reduce(self):
        self.chrom = ['chr' + c for c in [str(x) for x in range(1, 23)]]
        self.chrom.append('chrx')
        self.chrom.append('chry')
        self.export_reference_genome()
        self.export_reference_jellyfish_table()

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

    def export_jellyfish_table(self, channel, *args):
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
    # Simulation
    # ============================================================================================================================ #

    #def prepare(self):
    #    print('Simulation completed. Preparing statistics...')

    #def load_inputs(self):
    #    self.tracks_00 = {re.sub(r'\s+', '_', str(track).strip()).strip(): track for track in pybedtools.BedTool(os.path.join(self.get_current_job_directory(), 'absent.bed'))}
    #    self.tracks_10 = {re.sub(r'\s+', '_', str(track).strip()).strip(): track for track in pybedtools.BedTool(os.path.join(self.get_current_job_directory(), 'heterozygous.bed'))}
    #    self.tracks_11 = {re.sub(r'\s+', '_', str(track).strip()).strip(): track for track in pybedtools.BedTool(os.path.join(self.get_current_job_directory(), 'homozygous.bed'))}
    #    self.tracks = {}
    #    self.tracks.update(self.tracks_00)
    #    self.tracks.update(self.tracks_10)
    #    self.tracks.update(self.tracks_11)
    #    self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
    #    self.counts_provider = counttable.JellyfishCountsProvider(os.path.join(self.get_current_job_directory(), 'test.jf'))
    #    self.round_robin(self.tracks, filter_func = lambda track: track.end - track.start > 1000000) 

    #def transform(self, track, track_name):
    #    sv = self.get_sv_type()(track)
    #    c = config.Configuration()
    #    inner_kmers = sv.get_inner_kmers(self.reference_counts_provider.get_kmer_count, count = 1, n = 1000)
    #    if len(inner_kmers) == 0:
    #        return None
    #    path = os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track_name  + '.json')
    #    with open(path, 'w') as json_file:
    #        json.dump({'inner_kmers': {kmer: self.counts_provider.get_kmer_count(kmer) for kmer in inner_kmers} }, json_file, sort_keys = True, indent = 4)
    #    return path

    #def reduce(self):
    #    self.tracks = map_reduce.Job.reduce(self)
    #    self.kmers_00 = {}
    #    self.kmers_10 = {}
    #    self.kmers_11 = {}
    #    means_00 = []
    #    means_10 = []
    #    means_11 = []
    #    print(cyan(len(self.tracks)), 'tracks')
    #    for track in self.tracks:
    #        with open(self.tracks[track], 'r') as json_file:
    #            inner_kmers = json.load(json_file)['inner_kmers']
    #            #print(len(inner_kmers))
    #            if track in self.tracks_00:
    #                means_00.append(statistics.mean([inner_kmers[inner_kmer] for inner_kmer in inner_kmers]))
    #                for inner_kmer in inner_kmers:
    #                    self.kmers_00[inner_kmer] = inner_kmers[inner_kmer]
    #            if track in self.tracks_10:
    #                means_10.append(statistics.mean([inner_kmers[inner_kmer] for inner_kmer in inner_kmers]))
    #                for inner_kmer in inner_kmers:
    #                    self.kmers_10[inner_kmer] = inner_kmers[inner_kmer]
    #            if track in self.tracks_11:
    #                means_11.append(statistics.mean([inner_kmers[inner_kmer] for inner_kmer in inner_kmers]))
    #                for inner_kmer in inner_kmers:
    #                    self.kmers_11[inner_kmer] = inner_kmers[inner_kmer]
    #    print(green(len(self.tracks_00)), '00 tracks')
    #    print(green(len(self.tracks_10)), '10 tracks')
    #    print(green(len(self.tracks_11)), '11 tracks')
    #    print(blue(len(self.kmers_00)), '00 kmers')
    #    print(blue(len(self.kmers_10)), '10 kmers')
    #    print(blue(len(self.kmers_11)), '11 kmers')
    #    print(cyan(len(means_00)), '00 means')
    #    print(cyan(len(means_10)), '10 means')
    #    print(cyan(len(means_11)), '11 means')
    #    visualizer.histogram([self.kmers_00[kmer] for kmer in self.kmers_00], '00 kmers', self.get_current_job_directory(), x_label = 'number of kmers', y_label = 'coverage', step = 1)
    #    visualizer.histogram([self.kmers_10[kmer] for kmer in self.kmers_10], '10 kmers', self.get_current_job_directory(), x_label = 'number of kmers', y_label = 'coverage', step = 1)
    #    visualizer.histogram([self.kmers_11[kmer] for kmer in self.kmers_11], '11 kmers', self.get_current_job_directory(), x_label = 'number of kmers', y_label = 'coverage', step = 1)
    #    visualizer.histogram(means_00, 'mean coverage of 00 events', self.get_current_job_directory(), x_label = 'mean coverage', y_label = 'events', step = 1)
    #    visualizer.histogram(means_10, 'mean coverage of 10 events', self.get_current_job_directory(), x_label = 'mean coverage', y_label = 'events', step = 1)
    #    visualizer.histogram(means_11, 'mean coverage of 11 events', self.get_current_job_directory(), x_label = 'mean coverage', y_label = 'events', step = 1)

    #def plot_reference_kmer_profile(self):
    #    self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
    #    counts = []
    #    n = 0
    #    for kmer, count in self.reference_counts_provider.stream_kmers():
    #        print(kmer, count)
    #        counts.append(count)
    #        n += 1
    #        if n % 10000 == 1:
    #            print(green(n), 'kmers counted')
    #    counts = [counts[i] for i in random.sample(range(0, len(counts)), len(counts) / 10000)]
    #    visualizer.histogram(counts, 'reference_genome_kmer_profile', self.get_current_job_directory(), 'kmer coverage', 'number of kmers')
    #    exit()  

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
