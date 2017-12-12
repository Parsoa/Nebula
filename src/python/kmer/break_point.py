import io
import os
import re
import pwd
import sys
import copy
import json
import time
import argparse
import traceback

from multiprocessing import Process

from kmer import (
    bed,
    sets,
    config,
    commons,
    counttable,
    map_reduce,
    count_server,
)

from kmer.sv import StructuralVariation, Inversion, Deletion

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class BreakPoint(object):

    @staticmethod
    def to_json(break_point):
        return {
            'boundary': break_point.boundary,
            'kmers': break_point.kmers,
            'reference_kmers': break_point.reference_kmers
        }

    def __init__(self, boundary, begin, end, kmers, reference_kmers):
        self.name = '(' + str(begin) + ',' + str(end) + ')'
        self.boundary = boundary
        self.begin = begin
        self.end = end
        self.kmers = kmers
        self.reference_kmers = reference_kmers
        self.score = 0
        self.zygosity = None

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job for finding StructuralVariation breakpoints
# Algorithm: starts with a set of structural variation events and their approximate breakpoints and tries to refine them
# considers a radius of [-50, 50] around each end of the breakpoint, and for each pair of endpoints within that radius considers
# the area as the structural variations and applies it to the reference genome to generate a set of kmers. Discards those endpoints
# whose kmers do not all appear in the base genome the event was detected in. 
# Output: Reports the remaining boundary candidates with their list of associated kmers and the count of those kmers in the
# base genome.
# ============================================================================================================================ #
# ============================================================================================================================ #

class BreakPointJob(map_reduce.Job):

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def find_thread_count(self):
        pass

    def load_inputs(self):
        c = config.Configuration()
        bedtools = pybedtools.BedTool(c.bed_file)
        self.radius = 50
        # split variations into batches
        n = 0
        for track in bedtools:
            name = re.sub(r'\s+', '_', str(track).strip()).strip()
            # too large, skip
            if track.end - track.start > 1000000:
                print(colorama.Fore.RED + 'skipping ', name, ', too large')
                continue
            index = n % c.max_threads 
            if not index in self.batch:
                self.batch[index] = []
            self.batch[index].append(track)
            print(colorama.Fore.BLUE + 'assigned ', name, ' to ', index)
            n = n + 1
        self.num_threads = len(self.batch)
        print('running on ', self.num_threads, ' threads')

    def run_batch(self, batch):
        c = config.Configuration()
        sv_type = self.get_sv_type()
        output = {}
        for track in batch:
            name = re.sub(r'\s+', '_', str(track).strip()).strip()
            sv = sv_type(track = track, radius = self.radius)
            output[name] = self.transform(sv)
        self.output_batch(output)
        print(colorama.Fore.GREEN, 'process ', self.index, ' done')

    def transform(self, sv):
        c = config.Configuration()
        frontier = self.extract_boundary_kmers(sv)
        # whatever that is left in the frontier is a possible break point
        frontier = self.prune_boundary_candidates(frontier, sv)
        # now check the reference counts to find the best match
        results = {}
        results['candidates'] = len(frontier)
        for break_point in frontier :
            for kmer in break_point.reference_kmers:
                # counts for reference not available at this
                break_point.reference_kmers[kmer] = -1
            for kmer in break_point.kmers:
                break_point.kmers[kmer] = count_server.get_kmer_count(kmer, self.index, False)
            results[break_point.name] = BreakPoint.to_json(break_point)
            # save the number of boundary candidates
        return results

    # ============================================================================================================================ #
    # job-specific helpers
    # ============================================================================================================================ #

    def get_sv_type(self):
        c = config.Configuration()
        bed_file_name = c.bed_file.split('/')[-1]
        sv_type = bed_file_name.split('.')[-2]
        if sv_type == 'DEL':
            return Deletion
        if sv_type == 'INV':
            return Inversion
        return StructuralVariation

    def extract_boundary_kmers(self, sv):
        c = config.Configuration()
        frontier = {}
        for begin in range(-self.radius, self.radius + 1) :
            for end in range(-self.radius, self.radius + 1) :
                kmers, boundary = sv.get_signature_kmers(begin, end)
                if not kmers:
                    # skip this candidate
                    continue
                reference_kmers = sv.get_reference_signature_kmers(begin, end)
                #
                break_point = BreakPoint(boundary = boundary, begin = begin, end = end,\
                    kmers = kmers, reference_kmers = reference_kmers)
                frontier[break_point] = True
        return frontier

    # prunes a break points if not all its kmers appear in the counttable
    def prune_boundary_candidates(self, frontier, sv):
        c = config.Configuration()
        remove = {}
        for break_point in frontier:
            for kmer in break_point.kmers:
                count = count_server.get_kmer_count(kmer, self.index, False)
                if count == 0:
                    remove[break_point] = True
                    break
        for break_point in remove:
            frontier.pop(break_point, None)
        return frontier

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference")
    parser.add_argument("--bed")
    args = parser.parse_args()
    # 
    config.configure(reference_genome = args.reference, bed_file = args.bed)
    #
    break_point_job = BreakPointJob(job_name = 'break_point_', previous_job_name = '')
    break_point_job.execute()
