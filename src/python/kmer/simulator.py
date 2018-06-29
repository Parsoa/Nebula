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

    def execute(self):
        c = config.Configuration()
        self.check_cli_arguments(None)
        self.create_output_directories()
        self.simulate()
        exit()
        self.find_thread_count()
        self.prepare()
        self.load_inputs()
        if not self.resume_from_reduce:
            print('normal execution flow')
            self.distribute_workload()
            self.wait_for_children()
        else:
            print('resuming from reduce')
        output = self.reduce()
        self.plot(output)

    def simulate(self):
        c = config.Configuration()
        self.ref = extract_chromosome('chr1')
        if c.seed:
            random.seed(c.seed)
        if c.random:
            SVs = self.generate_random_intervals(ref, 1000)
        else:
            SVs = self.load_structural_variations()
        if c.heterozygous:
            pid = os.fork()
            # generate homozygous control channel using all events
            if pid == 0:
                self.export_bed(SVs, 'all')
                #control = self.generate_fasta(SVs)
                #self.export_fasta(control, 'control')
                #self.export_fastq('control')
                #self.export_jellyfish_table('control', 'control')
            # generate hetereozygous test channel using half the events
            else:
                n = len(SVs) / 2
                r = sorted(random.sample(xrange(len(SVs)), n))
                a = list(filter(lambda i: i not in r, range(len(SVs))))
                present = [ SVs[i] for i in r ]
                absent = [ SVs[i] for i in a ]
                self.export_bed(present, 'present')
                self.export_bed(absent, 'absent')
                strand_1, strand_2, homozygous_SVs, heterozygous_SVs = self.generate_fasta(present)
                self.export_bed(homozygous_SVs, 'homozygous')
                self.export_bed(heterozygous_SVs, 'heterozygous')
                self.export_fasta(strand_1, 'test_strand_1')
                self.export_fasta(strand_2, 'test_strand_2')
                # export a couple of fastq files for each strand, will this work?
                self.export_fastq('test_strand_1')
                self.export_fastq('test_strand_2')
                self.merge_fastq_files('test_strand_1')
                self.merge_fastq_files('test_strand_2')
                self.export_jellyfish_table('test', 'test_strand_1', 'test_strand_2')

    def generate_random_intervals(self, n):
        c = config.Configuration()
        random.seed(c.seed)
        print('generating', n, 'random events')
        l = len(self.ref)
        SVs = []
        offset = 100000
        for i in range(0, n):
            start = random.randint(offset, l - offset)
            end = start + random.randint(100, 2000)
            event = pybedtools.Interval(chrom = 'chr1', start = start, end = end)
            SVs.append(event)
        SVs = sorted(SVs, key = lambda x: x.start)
        SVs = self.filter_overlapping_intervals(SVs)
        return SVs
    
    def load_structural_variations(self):
        c = config.Configuration()
        bedtools = pybedtools.BedTool(c.bed_file)
        # split events into batches
        n = 0
        tracks = []
        for track in bedtools:
            if track.chrom != 'chr1':
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

    def generate_fasta(self, SVs):
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
        strand_1 = self.apply_events_to_ref(SVs)
        strand_2 = self.apply_events_to_ref(homozygous)
        return strand_1, strand_2, homozygous, heterozygous
    
    def export_fasta(self, sequence, name):
        print('Exporting FASTA file:', green(name + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), name + '.fa'), 'w') as fasta_file:
            fasta_file.write('>' + 'chr1' + '\n')
            l = len(sequence)
            num_lines = l / 50
            for i in range(0, num_lines):
                fasta_file.write(sequence[i * 50 : (i + 1) * 50] + '\n')
            if l % 50 != 0:
                fasta_file.write(sequence[num_lines * 50 :] + '\n')
        with open(os.path.join(self.get_current_job_directory(), name + '.lines.fa'), 'w') as fasta_lines_file:
            fasta_lines_file.write('>' + 'chr1' + '\n')
            fasta_lines_file.write(sequence)
    
    def apply_events_to_ref(self, SVs):
        seq = ''
        previous = 0
        for i, sv in enumerate(SVs):
            left = sv.start
            seq += self.ref[previous:left]
            if self.get_sv_type() == 'INV':
                seq += reverse_complement(self.ref[left:right])
            previous = sv.end
        return seq
    
    def export_fastq(self, name):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating FASTQ files:', green(name + '.1.fq'), green(name + '.2.fq'))
        num_reads = len(self.ref) * c.coverage / 100
        fasta = os.path.join(self.get_current_job_directory(), name + '.fa')
        fastq_1 = os.path.join(self.get_current_job_directory(), name + '.1.fq')
        fastq_2 = os.path.join(self.get_current_job_directory(), name + '.2.fq')
        command =  "wgsim -d400 -N{} -1100 -2100 {} {} {}".format(num_reads, fasta, fastq_1, fastq_2)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def merge_fastq_files(self, name):
        print('Merging fastq files', name)
        with open(os.path.join(self.get_current_job_directory(), name + '.1.fq'), 'a') as fastq_file:
                with open(os.path.join(self.get_current_job_directory(), name + '.2.fq'), 'r') as in_file:
                    for line in in_file:
                        fastq_file.write(line)
        os.rename(os.path.join(self.get_current_job_directory(), name + '.1.fq'), os.path.join(self.get_current_job_directory(), name + '.fq'))
        os.remove(os.path.join(self.get_current_job_directory(), name + '.2.fq'))

    #def export_bam(self, name):
        # generate sam file
        #print('Generating SAM file:', green(name + '.sam'))
        #sam_file = os.path.join(self.get_current_job_directory(), name + '.sam')
        #command = "bwa mem -M -t 20 -R \"@RG\\tID:1\\tPL:ILLUMINA\\tSM:cnv_1000_ref\" {} {} {} > {}".format(c.reference_genome, fastq_1, fastq_2, sam_file)
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        ## generate bam file
        #print('Generating unsorted BAM file:', green(name + '.unsorted.bam'))
        #unsorted_bam = os.path.join(self.get_current_job_directory(), name + '.unsorted.bam')
        #command = "samtools view -S -b {} > {}".format(sam_file, unsorted_bam)  
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        #print('Sorting ...')
        #bam = os.path.join(self.get_current_job_directory(), name)
        #command = "samtools sort {} -o {}".format(unsorted_bam, bam)
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        #print('Indexing ...')
        #bam_index = os.path.join(self.get_current_job_directory(), name + '.bai')
        #command = "samtools index {} {}".format(bam + '.bam', bam_index)
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        #print('Done!')
    
    def export_jellyfish_table(self, channel, *args):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating Jellyfish table')
        command = "jellyfish count -m 31 -s 1000000000 -t 24 --canonical "
        for name in args:
            command += ' ' + os.path.join(self.get_current_job_directory(), name + '.fq')
        command += ' -o ' + os.path.join(self.get_current_job_directory(), channel + '.jf')
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    # ============================================================================================================================ #
    # Simulation
    # ============================================================================================================================ #

    def prepare(self):
        print('Simulatio completed. Preparing statistics...')

    def load_inputs(self):
        self.tracks_00 = {re.sub(r'\s+', '_', str(track).strip()).strip(): track for track in pybedtools.BedTool(os.path.join(self.get_current_job_directory(), 'absent.bed'))}
        self.tracks_10 = {re.sub(r'\s+', '_', str(track).strip()).strip(): track for track in pybedtools.BedTool(os.path.join(self.get_current_job_directory(), 'heterozygous.bed'))}
        self.tracks_11 = {re.sub(r'\s+', '_', str(track).strip()).strip(): track for track in pybedtools.BedTool(os.path.join(self.get_current_job_directory(), 'homozygous.bed'))}
        self.tracks = {}
        self.tracks.update(self.tracks_00)
        self.tracks.update(self.tracks_10)
        self.tracks.update(self.tracks_11)
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish[0])
        self.counts_provider = counttable.JellyfishCountsProvider(os.path.join(self.get_current_job_directory(), 'test.jf'))
        self.round_robin(self.tracks, filter_func = lambda track: track.end - track.start > 1000000) 

    def transform(self, track, track_name):
        sv = self.get_sv_type()(track)
        c = config.Configuration()
        inner_kmers = sv.get_inner_kmers(self.reference_counts_provider.get_kmer_count, count = 1, n = 1000)
        if len(inner_kmers) == 0:
            return None
        path = os.path.join(self.get_current_job_directory(), 'inner_kmers_' + track_name  + '.json')
        with open(path, 'w') as json_file:
            json.dump({'inner_kmers': {kmer: self.counts_provider.get_kmer_count(kmer) for kmer in inner_kmers} }, json_file, sort_keys = True, indent = 4)
        return path

    def reduce(self):
        self.tracks = map_reduce.Job.reduce(self)
        self.kmers_00 = {}
        self.kmers_10 = {}
        self.kmers_11 = {}
        means_00 = []
        means_10 = []
        means_11 = []
        print(cyan(len(self.tracks)), 'tracks')
        for track in self.tracks:
            with open(self.tracks[track], 'r') as json_file:
                inner_kmers = json.load(json_file)['inner_kmers']
                #print(len(inner_kmers))
                if track in self.tracks_00:
                    means_00.append(statistics.mean([inner_kmers[inner_kmer] for inner_kmer in inner_kmers]))
                    for inner_kmer in inner_kmers:
                        self.kmers_00[inner_kmer] = inner_kmers[inner_kmer]
                if track in self.tracks_10:
                    means_10.append(statistics.mean([inner_kmers[inner_kmer] for inner_kmer in inner_kmers]))
                    for inner_kmer in inner_kmers:
                        self.kmers_10[inner_kmer] = inner_kmers[inner_kmer]
                if track in self.tracks_11:
                    means_11.append(statistics.mean([inner_kmers[inner_kmer] for inner_kmer in inner_kmers]))
                    for inner_kmer in inner_kmers:
                        self.kmers_11[inner_kmer] = inner_kmers[inner_kmer]
        print(green(len(self.tracks_00)), '00 tracks')
        print(green(len(self.tracks_10)), '10 tracks')
        print(green(len(self.tracks_11)), '11 tracks')
        print(blue(len(self.kmers_00)), '00 kmers')
        print(blue(len(self.kmers_10)), '10 kmers')
        print(blue(len(self.kmers_11)), '11 kmers')
        print(cyan(len(means_00)), '00 means')
        print(cyan(len(means_10)), '10 means')
        print(cyan(len(means_11)), '11 means')
        visualizer.histogram([self.kmers_00[kmer] for kmer in self.kmers_00], '00 kmers', self.get_current_job_directory(), x_label = 'number of kmers', y_label = 'coverage', step = 1)
        visualizer.histogram([self.kmers_10[kmer] for kmer in self.kmers_10], '10 kmers', self.get_current_job_directory(), x_label = 'number of kmers', y_label = 'coverage', step = 1)
        visualizer.histogram([self.kmers_11[kmer] for kmer in self.kmers_11], '11 kmers', self.get_current_job_directory(), x_label = 'number of kmers', y_label = 'coverage', step = 1)
        visualizer.histogram(means_00, 'mean coverage of 00 events', self.get_current_job_directory(), x_label = 'mean coverage', y_label = 'events', step = 1)
        visualizer.histogram(means_10, 'mean coverage of 10 events', self.get_current_job_directory(), x_label = 'mean coverage', y_label = 'events', step = 1)
        visualizer.histogram(means_11, 'mean coverage of 11 events', self.get_current_job_directory(), x_label = 'mean coverage', y_label = 'events', step = 1)

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
