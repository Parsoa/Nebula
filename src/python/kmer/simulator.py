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
        if c.diploid:
            print('Starting diploid simulation...')
        #
        self.extract_whole_genome()
        self.round_robin(self.chroms)

    def extract_whole_genome(self):
        self.chroms = extract_whole_genome()
        for chrom, seq in self.chroms.iteritems():
            print(chrom, 'with', len(seq), 'bases')
            with open(os.path.join(self.get_current_job_directory(), chrom + '.fa'), 'w') as chrom_file:
                chrom_file.write('>' + chrom + '\n')
                chrom_file.write(seq.strip())
                chrom_file.write('\n')

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
        self.sv_type = self.get_sv_type()
        print('Type of events:', cyan(self.sv_type))
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

    def get_sv_type(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        if bed_file_name.find('DEL') != -1:
            return 'DEL'
        if bed_file_name.find('INV') != -1:
            return 'INV'
        if bed_file_name.find('ALU') != -1:
            return 'ALU'
        return Deletion

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
                        track.zygosity = '0|0'
                    else:
                        self.present.append(track)
                        if r == 1:
                            self.heterozygous.append(track)
                            track.zygosity = '0|1'
                        else:
                            self.homozygous.append(track)
                            track.zygosity = '1|1'
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
        c = config.Configuration()
        strand_1 = self.apply_events_to_chromosome(chrom, self.present)
        strand_2 = self.apply_events_to_chromosome(chrom, self.homozygous)
        if c.diploid:
            self.export_diploid_chromosome_fasta(chrom, [strand_1, strand_2])
            exit()
            self.export_fastq(seq, chrom + '_diploid')
        else:
            self.export_chromosome_fasta(chrom, strand_1, chrom + '_strand_1')
            self.export_chromosome_fasta(chrom, strand_2, chrom + '_strand_2')
            exit()
            self.export_fastq(seq, chrom + '_strand_1')
            self.export_fastq(seq, chrom + '_strand_2')
        exit()

    def export_diploid_chromosome_fasta(self, chrom, strands):
        print('Exporting FASTA file:', green(chrom + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), chrom + '_diploid.fa'), 'w') as fasta_file:
            for i in range(0, 2):
                seq = strands[i]
                fasta_file.write('>' + chrom + '_' + str(i + 1) + '\n')
                fasta_file.write(seq.strip())
                fasta_file.write('\n')
                continue
                l = len(seq)
                n = 100
                num_lines = l / n
                for i in range(0, num_lines):
                    line = seq[i * n : (i + 1) * n].upper() + '\n'
                    fasta_file.write(line)

    def export_chromosome_fasta(self, chrom, seq, name):
        print('Exporting FASTA file:', green(name + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), name + '.fa'), 'w') as fasta_file:
            fasta_file.write('>' + chrom + '\n')
            fasta_file.write(seq.strip())
            fasta_file.write('\n')
            return
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
            if track.chrom.lower() != chrom.lower():
                print(yellow('skipping', track, 'on', chrom))
                continue
            print('applying', track, 'to', chrom)
            s = self.chroms[chrom][previous:track.begin]
            seq += s
            if self.sv_type == 'ALU':
                print('Adding sequence')
                seq += 'GCCGGGCGTGGTGGCTTACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGACGGGTGGATCACGAGGTCAGCAGATGGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATGCAAAAAAATTAGCCGGGTGTGGTGGTGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATGGCATGAACCTGGGAGGCGGAGCTTGCAGTGAGCCGAGATCATGTCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCGTCTCAAAAAAAAAAAGAAAAAAA'
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
        c = config.Configuration()
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
        command += ' -o ' + os.path.join(self.get_current_job_directory(), 'reference_' + str(c.ksize) + 'k.jf')
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def merge_fastq_files(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Mixing FASTQ files...')
        command = 'cat '
        if c.diploid:
            command += os.path.join(self.get_current_job_directory(), 'chr*.1.fq')
            command += ' > '
            command += os.path.join(self.get_current_job_directory(), 'test.1.fq')
            print(command)
            #output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)
            command = 'cat '
            command += os.path.join(self.get_current_job_directory(), 'chr*.2.fq')
            command += ' > '
            command += os.path.join(self.get_current_job_directory(), 'test.2.fq')
            print(command)
            #output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)
        else:
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

    def export_bam(name, ref):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        fasta = os.path.join(get_output_directory(), 'reference.fa')
        fastq = os.path.join(get_output_directory(), 'test.fq')
        # generate sam file
        print('Generating SAM file')
        sam_file = os.path.join(get_output_directory(), 'test.sam')
        command = "bwa mem -M -t 24 {} {} > {}".format(fasta, fastq, sam_file)
        subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        # generate bam file
        print('Generating unsorted BAM file')
        unsorted_bam = os.path.join(get_output_directory(), name + '.unsorted.bam')
        command = "samtools view -S -b {} > {}".format(sam_file, unsorted_bam)  
        subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        #print('Sorting ...')
        #bam = os.path.join(get_output_directory(), name)
        #command = "samtools sort {} -o {}".format(unsorted_bam, bam)
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        #print('Indexing ...')
        #bam_index = os.path.join(get_output_directory(), name + '.bai')
        #command = "samtools index {} {}".format(bam + '.bam', bam_index)
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        print('Done!')

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.resume_from_reduce)
