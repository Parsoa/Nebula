import io
import os
import re
import pwd
import sys
import copy
import json
import time
import traceback

from kmer import (
    bed,
    sets,
    config,
    commons,
    reference,
    counttable,
    count_server,
)

from kmer.sv import StructuralVariation

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class Genotype(object):

    def __init__(self, head, tail):
        self.head = head
        self.tail = tail
        self.kmers = bed.count_boundary_kmers(head, tail)
        self.score = 0
        self.kmer_list = list(self.kmers.keys())
        # print('distinc kmers in genotype: ', len(self.kmers))

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

radius = 50

# @commons.measure_time
def refine_variation_boundaries():
    c = config.Configuration()
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
        break
    # run each batch in a separate process
    children = []
    for index in batch:
        tracks = batch[index]
        pid = os.fork()
        if pid == 0:
            # forked process
            run_batch(tracks, index)
        else:
            # main process
            children.append(pid)
            print('spawned child ', pid)

def run_batch(tracks, index):
    output = {}
    for track in tracks:
        name = re.sub(r'\s+', '_', str(track).strip()).strip()
        print(colorama.Fore.GREEN + '========================================================')
        print(colorama.Fore.GREEN + 'track: ', name, '@', index)
        sv = StructuralVariation(track = track, radius = radius)
        output[name] = find_track_boundaries(sv, index)
    print('process ', index, ' done')
    # output manually, io redirection could get entangled with multiple client/servers
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/batch_' + str(index) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    exit()

@commons.measure_time
def find_track_boundaries(sv , index):
    c = config.Configuration()
    frontier = {}
    for begin in range(-radius, radius + 1) :
        for end in range(-radius, radius + 1) :
            head, tail = sv.get_interval_boundaries(begin, end)
            genotype = Genotype(head = head, tail = tail)
            key = str(begin) + '_' + str(end)
            frontier[key] = genotype
    # there will be c.ksize kmers at max
    for i in range(0, c.ksize) :
        print('i = ', i)
        print(len(frontier))
        remove = {}
        for key in frontier :
            genotype = frontier[key]
            kmers = []
            n = 0
            if 2 * i < len(genotype.kmer_list):
                kmers.append(genotype.kmer_list[2 * i])
                n = 1
                if (2 * i) + 1 < len(genotype.kmer_list):
                    kmers.append(genotype.kmer_list[2 * i + 1])
                    n = 2
            else :
                # due to repeats, it is possible that less than 2*ksize unique kmers appear
                continue
            genotype.score += score
            score = calc_similarity_score(kmers, index)
            if score != n :
                remove[key] = True
        for genotype in remove:
            print('removed: ', genotype)
            frontier.pop(genotype, None)
    # check whatever that is left in the frontier
    for key in frontier:
        print(frontier[key].head, '...', frontier[key].tail)
    # now check the reference counts to find the best match
    max = None
    # for begin in range(-radius, radius + 1) :
    #     for end in range(-radius, radius + 1) :
    #         try :
    #             interval = (begin, end)
    #             print(interval)
    #             print(colorama.Fore.GREEN + '--------------------------------------------------------')
    #             head, tail = sv.get_interval_boundaries(begin, end)
    #             score, boundary = calc_boundary_score(head, tail, index)
    #             # print(colorama.Fore.GREEN + 'score: ', score)
    #             max = (score, boundary, interval) if not max else\
    #                 (score, boundary, interval) if score > max[0] else max
    #         except Exception as e:
    #             print(e)
    #             traceback.print_exc()
    #             print(colorama.Fore.RED + 'boundary error, skipping')
    # print(colorama.Fore.GREEN + '########################################################')
    print('choice: ', max)
    return max

@commons.measure_time
def calc_boundary_score(head, tail, index):
    # print(colorama.Fore.GREEN + '--------------------------------------------------------')
    # print(colorama.Fore.GREEN + 'range: [', boundary[0], ', ', boundary[1], ']')
    # track = copy.deepcopy(t)
    # track.start = track.start + boundary[0]
    # track.end   = track.end   + boundary[1]
    # print(colorama.Fore.GREEN + 'range: [', track.start, ', ', track.end, ']')
    # reference_boundaries, variation_boundaries = bed.extract_track_boundaries(track)
    # we are not interested in the reference hereq
    kmers = bed.count_boundary_kmers(head, tail)
    # print(variation_kmers)
    score = len(calc_similarity_score(kmers, index))
    # score = 0
    return score, (head, tail)

def calc_similarity_score(kmers, index):
    result = {}
    for kmer in kmers:
        count = count_server.get_kmer_count(kmer, index)
        if count :
            # print(kmer, '{:04d}'.format(counttable.get_kmer_counts(kmer)[0]))
            # print(kmer, '{:04d}'.format(count))
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
