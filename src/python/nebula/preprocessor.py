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
import subprocess

from nebula import (
    bed,
    vcf,
    config,
    counter,
    junction,
    reduction,
    map_reduce,
    statistics,
    visualizer,
    programming
)

from nebula.debug import *
from nebula.kmers import *
from nebula.logger import *
from nebula.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class EventSetUnificationJob(map_reduce.Job):

    _name = 'EventSetUnificationJob'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = EventSetUnificationJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    SUPPORTED_SVTYPES = ['<DEL>', '<INS>', '<INV>']

    def execute(self):
        c = config.Configuration()
        self.create_output_directories()
        tracks = {}
        all_tracks = {}
        for path in c.bed:
            tracks[path] = {}
            tracks[path].update(vcf.load_tracks_from_file_as_dict(path, parse_header = True))
            for track in tracks[path]:
                tracks[path][track]['genotypes'] = {path: tracks[path][track]['genotype'].replace('|', '/')}
                if track not in all_tracks:
                    all_tracks.update({track: tracks[path][track]})
                else:
                    all_tracks[track]['genotypes'].update(tracks[path][track]['genotypes'])
            user_print('Loaded', len(tracks[path]), 'tracks from', path)
        user_print('Loaded', len(all_tracks), 'tracks.')
        # add genotypes for all samples to each track
        for track in all_tracks:
            for path in tracks:
                if path not in all_tracks[track]['genotypes']:
                    all_tracks[track]['genotypes'].update({path: '0/0'})
            assert len(all_tracks[track]['genotypes']) == 3
        all_tracks = bed.sort_tracks(all_tracks)
        # sort each svtype
        svtypes = {track.alt: True for track in all_tracks}
        #svtypes = {track.svtype: True for track in all_tracks}
        print(svtypes)
        merged_tracks = []
        for svtype in svtypes:
            if svtype in EventSetUnificationJob.SUPPORTED_SVTYPES:
                print('Processing ' + svtype + '..')
                _tracks = [track for track in all_tracks if track.alt == svtype]
                print('Merging', svtype + '.', len(_tracks), 'tracks.')
                _tracks = self.filter_overlapping_tracks(_tracks, svtype)
                print('Merged.', len(_tracks), 'tracks remaining.')
                merged_tracks += _tracks
        # merge all svtypes back together
        merged_tracks = bed.sort_tracks(merged_tracks)
        files = {}
        for path in c.bed:
            name = path.split('/')[-1]
            name = name[0: name.rfind('.')]
            #name = name + '.unified.all.vcf'
            files[path] = [open(os.path.join(self.get_current_job_directory(), name + '.unified.vcf'), 'w'), open(os.path.join(self.get_current_job_directory(), name + '.unified.repeat.vcf'), 'w'), open(os.path.join(self.get_current_job_directory(), name + '.unified.all.vcf'), 'w')]
        n = 0
        for track in merged_tracks:
            # for Illumina tracks
            track.svtype = track.alt[1:-1]
            genotypes = copy.deepcopy(track.genotypes)
            del track.genotypes
            for path in files:
                if n == 0:
                    files[path][0].write(track.header())
                    files[path][1].write(track.header())
                    files[path][2].write(track.header())
                if genotypes[path] != '0/0':
                    track['genotype'] = genotypes[path]
                    if track['IS_TRF'] == 'TR':
                        files[path][1].write(track.serialize())
                    else:
                        files[path][0].write(track.serialize())
                    files[path][2].write(track.serialize())
            n += 1
        exit()

    def calc_checksum(self, l):
        s = 0
        n = 0
        p = 0
        for track in l:
            if n % 2 == 0:
                s += int(track.begin)
            else:
                s -= int(track.begin)
            p = int(track.begin)
            n += 1
        print('Checksum:', s)

    def filter_overlapping_tracks(self, tracks, svtype):
        i = 0
        remove = []
        tracks = bed.sort_tracks(tracks)
        while i < len(tracks):
            for j in range(i + 1, len(tracks)):
                if tracks[j].chrom != tracks[i].chrom:
                    i = j
                    break
                if svtype == '<DEL>':
                    if tracks[j].begin <= int(tracks[i].end):
                        remove.append(j)
                        for path in tracks[j]['genotypes']:
                            if tracks[j]['genotypes'][path] != '0/0':
                                tracks[i]['genotypes'][path] = tracks[j]['genotypes'][path]
                        user_print_error('DEL Merged', str(tracks[j]), 'into', blue(str(tracks[i])))
                        continue
                # consider these the same track
                if svtype == '<INS>':
                    if tracks[j].begin - int(tracks[i].begin) < 100:
                        remove.append(j)
                        for path in tracks[j]['genotypes']:
                            if tracks[j]['genotypes'][path] != '0/0':
                                tracks[i]['genotypes'][path] = tracks[j]['genotypes'][path]
                        user_print_warning('INS Merged', str(tracks[j]), 'into', blue(str(tracks[i])))
                        continue
                i = j
            if j == len(tracks) - 1:
                break
        n = 0
        for index in sorted(remove):
            tracks.pop(index - n)
            n = n + 1
        return tracks

    def get_current_job_directory(self):
        return self.get_output_directory()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class TrackPreprocessorJob(map_reduce.Job):

    _name = 'TrackPreprocessorJob'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = TrackPreprocessorJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    SUPPORTED_SVTYPES = ['DEL', 'INS', 'INV', 'ALU', 'MEI']

    def execute(self):
        c = config.Configuration()
        tracks = {}
        for path in c.bed:
            tracks.update(bed.load_tracks_from_file_as_dict(path, parse_header = True))
        user_print('Loaded', len(tracks), 'tracks.')
        tracks = bed.sort_tracks(tracks)
        svtypes = {track.svtype: True for track in tracks}
        merged_tracks = []
        for svtype in svtypes:
            if svtype in TrackPreprocessorJob.SUPPORTED_SVTYPES:
                _tracks = [track for track in tracks if track.svtype == svtype]
                print('Merging', svtype + '.', len(_tracks), 'tracks.')
                _tracks = bed.filter_overlapping_tracks(_tracks, svtype)
                print('Merged.', len(_tracks), 'tracks remaining.')
                merged_tracks += _tracks
        tracks = bed.sort_tracks(merged_tracks)
        tracks = {track.id: track for track in tracks if track.svtype in TrackPreprocessorJob.SUPPORTED_SVTYPES}
        user_print('Removed overlapping tracks.', len(tracks), 'non-overlapping tracks.')
        return tracks

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class GcContentKmerSelectionJob(map_reduce.GenomeDependentJob):

    _name = 'ChromosomeGcContentEstimationJob'
    _category = 'preprocessing'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = ChromosomeGcContentEstimationJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def plot(self):
        values = open(os.path.join(self.get_current_job_directory(), 'depth_all.bed'), 'r').readlines()
        values = [float(line) for line in values]
        visualizer.scatter(list(range(0, 100)), values, 'GC_content_bp', self.get_current_job_directory(), 'GC', 'Coverage') 
        exit()

    def load_inputs(self):
        c = config.Configuration()
        self.gc = {}
        self.window_size = 500
        self.bin_size = self.window_size // 100
        for i in range(0, 100 + 1):
            self.gc[i] = {}
        self.load_reference_counts_provider() 
        self.chromosomes = extract_whole_genome()
        self.load_exon_regions()
        self.round_robin(self.chromosomes)

    def load_exon_regions(self):
        c = config.Configuration()
        exons = bed.load_tracks_from_file_as_dict(c.bed[0])
        self.exons = {}
        for chrom in self.chromosomes:
            self.exons[chrom] = {}
        for exon in exons:
            chrom = exons[exon].chrom.lower()
            if chrom in self.exons:
                self.exons[chrom][exon] = exons[exon]
        print(len(self.exons))

    #def transform(self, sequence, chrom):
    #    c = config.Configuration()
    #    t = time.time()
    #    n = 0
    #    for exon in self.exons[chrom]:
    #        track = self.exons[chrom][exon]
    #        ends = track.exonends.split(',')
    #        begins = track.exonstarts.split(',')
    #        for j, begin in enumerate(begins[:-1]):
    #            seq = sequence[int(begin): int(ends[j])]
    #            l = len(seq)
    #            i = self.window_size / 2
    #            while i < l - self.window_size / 2:
    #                kmer = canonicalize(seq[i - c.ksize / 2: i + c.ksize / 2])
    #                if self.reference_counts_provider.get_kmer_count(kmer) == 1:
    #                    contig = 'chrX' if 'x' in chrom else 'chrY' if 'y' in chrom else chrom + ':' + str(int(begin) + i - self.window_size / 2) + '-' + str(int(begin) + i + self.window_size / 2)
    #                    command = "samtools depth -r %s %s | awk 'BEGIN {cnt = 0} {sum+=$3; cnt++} END {OFS = \"\t\"; print sum / cnt}'" % (contig, c.bam)
    #                    coverage = subprocess.Popen(command, stdout = subprocess.PIPE, shell = True).communicate()[0].strip()
    #                    try:
    #                        coverage = float(coverage)
    #                    except:
    #                        system_print_error(command)
    #                        i += 1
    #                        continue
    #                    gc = calculate_gc_content(seq[i - self.window_size / 2: i + self.window_size / 2]) / self.bin_size
    #                    self.gc[gc][kmer] = {'contig': contig, 'gc': gc, 'coverage': coverage}
    #                    i += self.window_size / 2 + c.ksize / 2
    #                else:
    #                    i += 1
    #                    continue
    #        n += 1
    #        if n % 100 == 0:
    #            s = time.time()
    #            p = n / float(len(self.exons[chrom]))
    #            e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
    #            print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
    #    return None

    # No exon, unique kmers from all across the genome
    def transform(self, sequence, chrom):
        c = config.Configuration()
        telomere_length = 10000
        sequence = sequence[telomere_length: -telomere_length]
        t = time.time()
        l = len(sequence)
        i = self.window_size // 2
        while i < l - self.window_size // 2:
            kmer = canonicalize(sequence[i - c.ksize // 2: i + c.ksize // 2])
            if self.reference_counts_provider.get_kmer_count(kmer) == 1:
                gc = calculate_gc_content(sequence[i - self.window_size // 2: i + self.window_size // 2]) // self.bin_size
                contig = chrom + ':' + str(telomere_length + i - self.window_size // 2) + '-' + str(telomere_length + i + self.window_size // 2)
                self.gc[gc][kmer] = {'contig': contig, 'gc': gc, 'coverage': gc}
                i += self.window_size // 2 + c.ksize // 2
            else:
                i += 1
                continue
            if i % 10000 == 0:
                s = time.time()
                p = i / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))
        return None

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.gc, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def reduce(self):
        c = config.Configuration()
        self.kmers = {}
        random.seed(c.seed)
        gc_files = []
        #for gc in range(0, 100 + 1):
        #    gc_files.append(open(os.path.join(self.get_current_job_directory(), 'gc_' + str(gc) + '.bed'), 'w'))
        for batch in self.load_output():
            for gc in range(0, 100 + 1):
                for kmer in batch[str(gc)]:
                    self.gc[gc][kmer] = batch[str(gc)][kmer]
                    #gc_files[gc].write(batch[str(gc)][kmer]['contig'] + '\n')
        for gc in self.gc:
            print(len(self.gc[gc]), 'kmers with GC content', gc)
            if len(self.gc[gc]) > 10000:
                k = list(self.gc[gc].keys())
                for kmer in [k[i] for i in sorted(random.sample(xrange(len(k)), 10000))]:
                    self.kmers[kmer] = self.gc[gc][kmer]
            else:
                for kmer in self.gc[gc]:
                    self.kmers[kmer] = self.gc[gc][kmer]
        print('counting', len(self.kmers), 'kmers')
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, indent = 4)
        return self.gc

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class MixKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'MixKmersJob'
    _category = 'preprocessing'
    _previous_job = None 
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = MixKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.load_kmers()
        self.merge_kmers()
        self.export_tracks()
        exit()

    def load_kmers(self):
        c = config.Configuration()
        self.load_inner_kmers()
        self.load_junction_kmers()
        self.load_gc_content_kmers()
        self.load_depth_of_coverage_kmers()
        system_print_high('Counting', green(len(self.inner_kmers)), 'inner kmers')
        system_print_high('Counting', green(len(self.junction_kmers)), 'junction kmers')
        system_print_high('Counting', green(len(self.depth_kmers)), 'depth-of-coverage kmers')
        system_print_high('Counting', green(len(self.gc_kmers)), 'GC-content kmers')

    def merge_kmers(self):
        kmers = {}
        kmers['gc_kmers'] = self.gc_kmers
        kmers['depth_kmers'] = self.depth_kmers
        kmers['inner_kmers'] = self.inner_kmers
        kmers['junction_kmers'] = self.junction_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

    def load_inner_kmers(self):
        self.inner_kmers = {}
        job = reduction.FilterInnerKmersJob()
        with open(os.path.join(job.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
            self.inner_kmers.update(json.load(json_file))

    def load_junction_kmers(self):
        c = config.Configuration()
        job = junction.FilterJunctionKmersJob()
        self.junction_kmers = {}
        with open(os.path.join(job.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                kmer = canonicalize(kmer)
                if kmer in self.inner_kmers:
                    self.inner_kmers.pop(kmer, None)
                else:
                    self.junction_kmers[kmer] = kmers[kmer]

    def load_depth_of_coverage_kmers(self):
        n = 100000
        self.depth_kmers = {}
        self.load_reference_counts_provider() 
        for kmer, count in self.reference_counts_provider.stream_kmers():
            if count == 1 and kmer.find('N') == -1:
                canon = canonicalize(kmer)
                if not canon in self.inner_kmers and not canon in self.junction_kmers:
                    self.depth_kmers[canon] = {'loci': {}, 'count': 0}
                    n -= 1
                    if n == 0:
                        break
        self.unload_reference_counts_provider()

    def load_gc_content_kmers(self):
        c = config.Configuration()
        job = GcContentKmerSelectionJob()
        gc_kmers = json.load(open(c.gckmers)) if c.gckmers else job.execute()
        self.gc_kmers = {}
        for kmer in gc_kmers:
            if kmer not in self.inner_kmers and kmer not in self.junction_kmers:
                self.gc_kmers[kmer] = gc_kmers[kmer]
                self.gc_kmers[kmer].update({'loci': {}, 'count': 0})

    def export_tracks(self):
        c = config.Configuration()
        self.tracks = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.junction_kmers:
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        system_print_high('Kmers exported for', len(self.tracks), 'tracks.')
        with open(os.path.join(self.get_current_job_directory(), 'tracks.bed'), 'w') as bed_file:
            for track in self.tracks:
                bed_file.write(c.tracks[track].serialize())
