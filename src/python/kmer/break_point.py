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
        self.kmer_list = list(self.kmers.keys())
        self.score = 0
        self.zygosity = None

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

radius = 50

@commons.measure_time
def execute():
    c = config.Configuration()
    bedtools = pybedtools.BedTool(c.bed_file)
    # split variations into batches
    n = 0
    batch = {}
    for track in bedtools:
        name = re.sub(r'\s+', '_', str(track).strip()).strip()
        if track.end - track.start > 1000000:
            print(colorama.Fore.RED, 'skipping ', name, ', too large')
            # too large skip
            continue
        index = n % c.num_threads 
        if not index in batch :
            batch[index] = []
        batch[index].append(track)
        print(colorama.Fore.BLUE, 'assigned ', name, ' to ', index)
        n = n + 1
    # run each batch in a separate process
    children = {}
    for index in batch:
        tracks = batch[index]
        pid = os.fork()
        if pid == 0:
            # forked process
            run_batch(tracks, index)
        else:
            # main process
            children[pid] = True
            print(colorama.Fore.BLUE, 'spawned child ', pid)
    while True:
        (pid, e) = os.wait()
        children.pop(pid, None)
        print(colorama.Fore.RED, 'pid ', pid, 'finished')
        if len(children) == 0:
            break
    print(colorama.Fore.BLUE, 'all children done, merging output', pid)
    # merge_outputs()

def merge_outputs():
    c = config.Configuration()
    output = {}
    for i in range(0, c.num_threads):
        with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_' + str(i) + '.json')), 'r') as json_file:
            batch = json.load(json_file)
            output.update(batch)
    bed_file_name = c.bed_file.split('/')[-1]
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/boundaries_' + bed_file_name + '_' + str(c.ksize) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    
def clean_up():
    c = config.Configuration()
    # delete intermediate results
    for i in range(0, c.num_threads):
        # this might fail if there were less samples than threads
        try:
            os.remove(os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_' + str(i) + '.json')))
        except Exception as e:
            continue
    
def run_batch(tracks, index):
    c = config.Configuration()
    sv_class = StructuralVariation.get_sv_class()
    output = {}
    for track in tracks:
        name = re.sub(r'\s+', '_', str(track).strip()).strip()
        print(colorama.Fore.GREEN + '========================================================')
        print(colorama.Fore.GREEN + 'track: ', name, '@', index)
        sv = sv_class(track = track, radius = radius)
        output[name] = find_track_boundaries(sv, index)
    print(colorama.Fore.GREEN, 'process ', index, ' done')
    # output manually, io redirection could get entangled with multiple client/servers
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/batch_' + str(index) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=False, indent=4, separators=(',', ': '))
    exit()

# @commons.measure_time
def find_track_boundaries(sv , index):
    c = config.Configuration()
    frontier = extract_boundary_kmers(sv, index)
    # whatever that is left in the frontier is a possible break point
    frontier = prune_boundary_candidates(frontier, sv, index)
    # now check the reference counts to find the best match
    results = {}
    results['candidates'] = len(frontier)
    for break_point in frontier :
        for kmer in break_point.reference_kmers:
            break_point.reference_kmers[kmer] = count_server.get_kmer_count(kmer, index, False)
        for kmer in break_point.kmers:
            break_point.kmers[kmer] = count_server.get_kmer_count(kmer, index, False)
        results[break_point.name] = BreakPoint.to_json(break_point)
        # save the number of boundary candidates
    return results

def extract_boundary_kmers(sv, index):
    c = config.Configuration()
    frontier = {}
    for begin in range(-radius, radius + 1) :
        for end in range(-radius, radius + 1) :
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

def prune_boundary_candidates(frontier, sv, index):
    c = config.Configuration()
    # there will be c.ksize kmers at max
    for i in range(0, c.ksize) :
        # print('i = ', i)
        # print(len(frontier))
        remove = {}
        for break_point in frontier :
            kmers = []
            # due to repeats, it is possible that less than 2*ksize unique kmers appear
            n = 0
            if 2 * i < len(break_point.kmer_list):
                kmers.append(break_point.kmer_list[2 * i])
                n = 1
                if (2 * i) + 1 < len(break_point.kmer_list):
                    kmers.append(break_point.kmer_list[2 * i + 1])
                    n = 2
            else :
                continue
            score = calc_similarity_score(kmers, index)
            break_point.score += score
            if score != n :
                remove[break_point] = True
        for break_point in remove:
            # print('removed: ', break_point.name)
            frontier.pop(break_point, None)
    return frontier

def calc_similarity_score(kmers, index):
    result = {}
    for kmer in kmers:
        count = count_server.get_kmer_count(kmer, index, False)
        if count:
            result[kmer] = count
    return len(result)

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference")
    parser.add_argument("--type")
    parser.add_argument("--bed")
    args = parser.parse_args()
    # 
    config.configure(reference_genome = args.reference, bed_file = args.bed,\
        variation_type = args.type)
    #
    execute()
