import io
import os
import re
import pwd
import sys
import copy
import json
import time
import traceback

from multiprocessing import Process

from kmer import (
    bed,
    sets,
    config,
    commons,
    reference,
    counttable,
    count_server,
)

from kmer.sv import StructuralVariation, Inversion, Deletion

import khmer
import colorama
import pybedtools

c = config.Configuration()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class BreakPoint(object):

    @staticmethod
    def to_json(break_point):
        return {
            '0 head': break_point.head,
            '1 tail': break_point.tail,
            '2 kmers': break_point.kmers,
            '3 reference_kmers': break_point.reference_kmers
        }

    def __init__(self, head, tail, begin, end, kmers, reference_kmers):
        self.name = '(' + str(begin) + ',' + str(end) + ')'
        self.begin = begin
        self.end = end
        self.head = head
        self.tail = tail
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
def refine_variation_boundaries():
    bedtools = pybedtools.BedTool(c.bed_file)
    # split variations into batches
    n = 0
    batch = {}
    for track in bedtools:
        name = re.sub(r'\s+', '_', str(track).strip()).strip()
        index = n % c.num_threads 
        if not index in batch :
            batch[index] = []
        batch[index].append(track)
        print('assigned ', name, ' to ', index)
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
            print('spawned child ', pid)
    while True:
        (pid, e) = os.wait()
        children.pop(pid, None)
        print(colorama.Fore.RED, 'pid ', pid, 'finished')
        if len(children) == 0:
            break
    print('all children done, merging output', pid)
    merge_outputs()

def merge_outputs():
    output = {}
    for i in range(0, 12):
        with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_' + str(i) + '.json')), 'r') as json_file:
            batch = json.load(json_file)
            output.update(batch)
    bed_file_name = c.bed_file.split('/')[-1]
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/boundaries_' + bed_file_name + '_' + str(c.ksize) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    # delete intermediate results
    for i in range(0, 12):
        os.remove(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_' + str(i) + '.json')))
    
def run_batch(tracks, index):
    output = {}
    for track in tracks:
        name = re.sub(r'\s+', '_', str(track).strip()).strip()
        print(colorama.Fore.GREEN + '========================================================')
        print(colorama.Fore.GREEN + 'track: ', name, '@', index)
        sv = Inversion(track = track, radius = radius)
        output[name] = find_track_boundaries(sv, index)
    print(colorama.Fore.GREEN, 'process ', index, ' done')
    # output manually, io redirection could get entangled with multiple client/servers
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/batch_' + str(index) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    exit()

# @commons.measure_time
def find_track_boundaries(sv , index):
    frontier = extract_boundary_kmers()
    # whatever that is left in the frontier is a possible break point
    frontier = prune_boundary_candidates(frontier)
    # now check the reference counts to find the best match
    results = {}
    results['candidates'] = len(frontier)
    for break_point in frontier :
        for kmer in break_point.reference_kmers:
            break_point.reference_kmers[kmer] = count_server.get_kmer_count(kmer, index)
        for kmer in break_point.kmers:
            break_point.kmers[kmer] = count_server.get_kmer_count(kmer, index)
        results[break_point.name] = BreakPoint.to_json(break_point)
        # save the number of boundary candidates
    return results

def extract_boundary_kmers():
    frontier = {}
    for begin in range(-radius, radius + 1) :
        for end in range(-radius, radius + 1) :
            head, tail = sv.get_signature_kmers(begin, end, True)
            if not head:
                # skip this candidate
                continue
            #
            kmers = count_server.count_kmers_exact_list(head, tail)
            #
            ref_head, ref_tail = sv.get_reference_signature_kmers(begin, end)
            reference_kmers = count_server.count_kmers_exact_list(ref_head, ref_tail)
            #
            break_point = BreakPoint(head = head, tail = tail, begin = begin, end = end,\
                kmers = kmers, reference_kmers = reference_kmers)
            frontier[break_point] = True
    return frontier

def prune_boundary_candidates(frontier):
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

def prune_kmers(kmers):
    remove = {}
    for kmer in kmers:
        count = count_server.get_kmer_count(kmer, index, True)
        if count:
            remove.append(kmer)
    for kmer in remove:
        kmers.pop(kmer, None)
    return kmer

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

def execute():
    refine_variation_boundaries()

if __name__ == '__main__':
    config.configure()
    execute()
