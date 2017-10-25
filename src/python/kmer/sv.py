import io
import os
import pwd
import sys
import copy
import json
import time

from . import (
    bed,
    sets,
    config,
    commons,
    reference,
    counttable,
    count_server,
)

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

radius = 50

@commons.measure_time
def refine_variation_boundaries():
    c = config.Configuration()
    bedtools = pybedtools.BedTool(c.bed_file)
    # split variations into batches
    n = 0
    batch = {}
    for track in bedtools:
        if not n % c.server_count in jobs :
            batch[n % c.server_count] = []
        batch[n % c.server_count].append(track)
        n = n + 1
    # run each batch in a separate process
    for index in batch:
        tracks = batch[index]
        # create a child process
        pid = os.fork()
        if pid == 0:
            # child
            output = {}
            for track in tracks:
                output[str(track).strip()] = find_track_boundaries(track)
            # output manually, io redirection could get entangled with multiple client/servers
            with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_' + str(index) + '.json')), 'w') as json_file:
                json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))
            exit()

@commons.measure_time
def find_track_boundaries(track):
    print(colorama.Fore.GREEN + '========================================================')
    print(colorama.Fore.GREEN + 'track: ', str(track).strip())
    #
    max = None
    for start in range(-radius, radius + 1) :
        for end in range(-radius, radius + 1) :
            try :
                interval = (start, end)
                score, boundary = calc_boundary_score(track, interval)
                print(colorama.Fore.GREEN + 'score: ', score)
                max = (score, boundary, interval) if not max else\
                    (score, boundary, interval) if score > max[0] else max
            except Exception as e:
                print(e)
                print(colorama.Fore.RED + 'boundary error, skipping')
    print(colorama.Fore.GREEN + '########################################################')
    print('choice: ', max)
    return max

# @commons.measure_time
def calc_boundary_score(t, boundary):
    print(colorama.Fore.GREEN + '--------------------------------------------------------')
    print(colorama.Fore.GREEN + 'range: [', boundary[0], ', ', boundary[1], ']')
    track = copy.deepcopy(t)
    track.start = track.start + boundary[0]
    track.end   = track.end   + boundary[1]
    # print(colorama.Fore.GREEN + 'range: [', track.start, ', ', track.end, ']')
    reference_boundaries, variation_boundaries = bed.extract_track_boundaries(track)
    # we are not interested in the reference here
    variation_kmers = bed.count_boundary_kmers(variation_boundaries)
    # 
    # print(variation_kmers)
    score = len(calc_similarity_score(variation_kmers))
    # score = 0
    return (score, variation_boundaries)

def calc_similarity_score(kmers):
    result = {}
    for kmer in kmers:
        count = count_server.get_kmer_count(kmer)
        if count :
            # print(kmer, '{:04d}'.format(counttable.get_kmer_counts(kmer)[0]))
            print(kmer, '{:04d}'.format(count))
            result[kmer] = count
    return result

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

def execute():
    refine_variation_boundaries()

if __name__ == '__main__':
    config.configure()
    execute()
