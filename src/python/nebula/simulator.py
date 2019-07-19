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

from nebula import (
    bed,
    counttable,
    map_reduce,
    statistics,
    visualizer,
)

from nebula.kmers import *
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# Required arguments:
# --seed: random seed to use
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
        random.seed(c.seed)
        self.tracks = self.load_structural_variations()
        self.assign_event_zygosities()
        self.export_bed(self.tracks, 'all')
        self.export_bed(self.absent, 'absent')
        self.export_bed(self.present, 'present')
        self.export_bed(self.homozygous, 'homozygous')
        self.export_bed(self.heterozygous, 'heterozygous')
        self.extract_chromosomes()
        self.round_robin(self.chroms)

    def extract_chromosomes(self):
        self.chroms = extract_whole_genome()
        for chrom, seq in self.chroms.iteritems():
            print(chrom, 'with', len(seq), 'bases')
            with open(os.path.join(self.get_current_job_directory(), chrom + '.fa'), 'w') as chrom_file:
                chrom_file.write('>' + chrom + '\n')
                chrom_file.write(seq.strip())
                chrom_file.write('\n')

    def load_structural_variations(self):
        c = config.Configuration()
        tracks = [c.tracks[track] for track in c.tracks]
        return self.filter_overlapping_tracks(\
                    sorted(sorted(tracks, key = lambda x: x.begin), key = lambda y: y.chrom)\
                )

    def filter_overlapping_tracks(self, tracks):
        remove = []
        i = 0
        while i < len(tracks):
            for j in range(i + 1, len(tracks)):
                # j is contained inside i
                if tracks[j].chrom != tracks[i].chrom:
                    i = j
                    break
                if tracks[j].begin <= tracks[i].end:
                    remove.append(j)
                    print(red(str(tracks[j])), 'overlaps', blue(str(tracks[i])))
                    continue
                else:
                    i = j
                    break
            if i == len(tracks) - 1:
                break
        n = 0
        for index in sorted(remove):
            tracks.pop(index - n)
            n = n + 1
        return tracks

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
                if hasattr(track, 'genotype'):
                    if track.genotype == '10' or track.genotype == '01':
                        self.heterozygous.append(track)
                        self.present.append(track)
                    elif track.genotype == '11':
                        self.homozygous.append(track)
                        self.present.append(track)
                    else:
                        self.absent.append(track)
                else:
                    r = random.randint(0, 2)
                    if r == 2:
                        self.absent.append(track)
                        track.add_field('genotype', '00')
                    else:
                        self.present.append(track)
                        if r == 1:
                            self.heterozygous.append(track)
                            track.add_field('genotype', '10')
                        else:
                            self.homozygous.append(track)
                            track.add_field('genotype', '11')
        print(len(self.tracks), 'non-overlapping tracks')
        print(len(self.absent), '0|0')
        print(len(self.present), '1|1 or 1|0')
        print(len(self.homozygous), '1|1')
        print(len(self.heterozygous), '1|0')

    def export_bed(self, tracks, name):
        print('Exporting BED file:', green(name + '.bed'))
        with open(os.path.join(self.get_current_job_directory(), name + '.bed'), 'w') as bed_file:
            for index, track in enumerate(tracks):
                if index == 0:
                    bed_file.write(track.header())
                bed_file.write(track.serialize()) 

    def transform(self, seq, chrom):
        print('simulating', chrom)
        c = config.Configuration()
        strand_1 = self.apply_events_to_chromosome(chrom, self.present)
        strand_2 = self.apply_events_to_chromosome(chrom, self.homozygous)
        if c.diploid:
            self.export_diploid_chromosome_fasta(chrom, [strand_1, strand_2])
            self.export_fastq(seq, chrom + '_diploid')
        else:
            self.export_chromosome_fasta(chrom, strand_1, chrom + '_strand_1')
            self.export_chromosome_fasta(chrom, strand_2, chrom + '_strand_2')
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
            if track.svtype == 'ALU':
                seq += 'GCCGGGCGTGGTGGCTTACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGACGGGTGGATCACGAGGTCAGCAGATGGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATGCAAAAAAATTAGCCGGGTGTGGTGGTGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATGGCATGAACCTGGGAGGCGGAGCTTGCAGTGAGCCGAGATCATGTCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCGTCTCAAAAAAAAAAAGAAAAAAA'
            if track.svtype == 'INS':
                seq += track.seq.upper()
            if track.svtype == 'INV':
                seq += reverse_complement(self.chroms[chrom][track.begin: track.end])
            previous = track.end
        seq += self.chroms[chrom][previous:]
        return seq

    def export_fastq(self, seq, name):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating FASTQ files:', green(name + '.1.fq'), green(name + '.2.fq'))
        num_reads = len(seq) * c.coverage / 100
        fasta = os.path.join(self.get_current_job_directory(), name + '.fa')
        fastq_1 = os.path.join(self.get_current_job_directory(), name + '.1.fq')
        fastq_2 = os.path.join(self.get_current_job_directory(), name + '.2.fq')
        command =  "wgsim -d400 -N{} -1100 -2100 -r{} -R0 -e{} {} {} {}".format(num_reads, c.mutation_rate, c.sequencing_error_rate, fasta, fastq_1, fastq_2)
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def reduce(self):
        c = config.Configuration()
        self.export_reference_genome()
        #self.export_reference_jellyfish_table()
        #self.merge_fastq_files()

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
        command += os.path.join(self.get_current_job_directory(), 'chr*.1.fq')
        command += ' > '
        command += os.path.join(self.get_current_job_directory(), 'strand.1.fq')
        print(command)
        #output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)
        command = 'cat '
        command += os.path.join(self.get_current_job_directory(), 'chr*.2.fq')
        command += ' > '
        command += os.path.join(self.get_current_job_directory(), 'strand.2.fq')
        print(command)

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
