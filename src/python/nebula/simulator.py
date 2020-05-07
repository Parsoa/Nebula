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
)

from nebula.kmers import *
from nebula.logger import *
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
        for chrom, seq in sorted(self.chroms.items(), key = lambda x: x[0]):
            print(chrom, 'with', len(seq), 'bases')
            with open(os.path.join(self.get_current_job_directory(), chrom + '.fa'), 'w') as chrom_file:
                chrom_file.write('>' + chrom + '\n')
                chrom_file.write(seq.strip())
                chrom_file.write('\n')

    def load_structural_variations(self):
        c = config.Configuration()
        return self.filter_overlapping_tracks(c.tracks)

    def filter_overlapping_tracks(self, tracks):
        i = 0
        remove = []
        tracks = bed.sort_tracks(tracks)
        while i < len(tracks):
            for j in range(i + 1, len(tracks)):
                if tracks[j].chrom != tracks[i].chrom:
                    i = j
                    break
                if tracks[j].begin <= tracks[i].end:
                    remove.append(j)
                    user_print_warning(str(tracks[j]), 'overlaps', blue(str(tracks[i])))
                    continue
                if tracks[j].begin - tracks[i].end < 1000:
                    remove.append(j)
                    user_print_warning(str(tracks[j]), 'is too close to', blue(str(tracks[i])))
                    continue
                i = j
            if j == len(tracks) - 1:
                break
        n = 0
        for index in sorted(remove):
            tracks.pop(index - n)
            n = n + 1
        return tracks
    
    def assign_event_zygosities(self):
        c = config.Configuration()
        self.absent = []
        self.present = []
        self.homozygous = []
        self.heterozygous = []
        for track in self.tracks:
            r = random.randint(0, 2 if c.random else 1)
            if r == 2:
                self.absent.append(track)
                track['genotype'] = '0/0'
            else:
                chrom = track.chrom
                if chrom == 'chrx':
                    self.present.append(track)
                    if c.gender == 'Male':
                        self.heterozygous.append(track)
                        track['genotype'] = '1/0'
                    else:
                        if r == 1:
                            self.heterozygous.append(track)
                            track['genotype'] = '1/0'
                        else:
                            self.homozygous.append(track)
                            track['genotype'] = '1/1'
                elif chrom == 'chry':
                    if c.gender == 'Female':
                        self.present.append(track)
                        self.heterozygous.append(track)
                        track['genotype'] = '1/0'
                    else:
                        self.absent.append(track)
                        track['genotype'] = '0/0'
                else:
                    self.present.append(track)
                    if r == 1:
                        self.heterozygous.append(track)
                        track['genotype'] = '1/0'
                    else:
                        self.homozygous.append(track)
                        track['genotype'] = '1/1'
        print(len(self.tracks), 'non-overlapping tracks')
        print(len(self.absent), '0/0')
        print(len(self.present), '1/1 or 1/0')
        print(len(self.homozygous), '1/1')
        print(len(self.heterozygous), '1/0')

    def export_bed(self, tracks, name):
        print('Exporting BED file:', green(name + '.bed'))
        with open(os.path.join(self.get_current_job_directory(), name + '.bed'), 'w') as bed_file:
            bed_file.write(self.present[0].header())
            for index, track in enumerate(tracks):
                bed_file.write(track.serialize()) 

    def transform(self, seq, chrom):
        print('Simulating', chrom)
        c = config.Configuration()
        if chrom.lower() != 'chrx' and chrom.lower() != 'chry':
            strand_1 = self.apply_events_to_chromosome(chrom, self.present)
            strand_2 = self.apply_events_to_chromosome(chrom, self.homozygous)
        elif chrom.lower() == 'chrx':
            if c.gender == 'Male':
                strand_1 = self.apply_events_to_chromosome('chrX', self.present)
                strand_2 = self.apply_events_to_chromosome('chrY', self.homozygous)
            else:
                strand_1 = self.apply_events_to_chromosome('chrX', self.present)
                strand_2 = self.apply_events_to_chromosome('chrX', self.homozygous)
        else:
            # chrY
            return None
        self.export_diploid_chromosome_fasta(chrom, [strand_1, strand_2])
        self.export_fastq(seq, chrom + '_diploid')
        return None

    def apply_events_to_chromosome(self, chrom, tracks):
        seq = ''
        previous = 0
        for track in tracks:
            if track.chrom.lower() != chrom.lower():
                #print(yellow('skipping', track, 'on', chrom))
                continue
            #print('applying', track, 'to', chrom)
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

    def export_diploid_chromosome_fasta(self, chrom, strands):
        print('Exporting FASTA file:', green(chrom + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), chrom + '_diploid.fa'), 'w') as fasta_file:
            for i in range(0, 2):
                fasta_file.write('>' + chrom + '_' + str(i + 1) + '\n')
                fasta_file.write(strands[i])
                fasta_file.write('\n')
                continue

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

    def output_batch(self, batch):
        exit()

    def reduce(self):
        c = config.Configuration()
        self.export_reference_genome()

# ============================================================================================================================ #
# ============================================================================================================================ #
# Required arguments
# --seed: random seed to use
# ============================================================================================================================ #
# ============================================================================================================================ #

class KmerCountSimulator(map_reduce.Job):

    _name = 'KmerCountSimulator'
    _category = 'simulation'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = KmerCountSimulator(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # Simulation
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        path = c.kmers[0]
        with open(path, 'r') as json_file:
            self.kmers = json.load(json_file)
        for kmer_type in ['depth_kmers', 'inner_kmers', 'junction_kmers']:
            for kmer in kmers[kmer_type]:
                kmers[kmer_type][kmer] = None
