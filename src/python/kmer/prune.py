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
    break_point,
    count_server,
)

from kmer.sv import StructuralVariation, Inversion, Deletion

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

def execute():
    c = config.Configuration()
    # split cnadidates into batches
    batch = {}
    for index in range(0, c.num_threads):
        path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_' + str(index) + '.json'))
        if os.path.isfile(path):
            with open(path, 'r') as json_file:
                batch[index] = json.load(json_file)
    # run each batch
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
    print('all forks done, merging output ...', pid)
    merge_outputs()

def merge_outputs():
    c = config.Configuration()
    output = {}
    for i in range(0, c.num_threads):
        # might fail because there weren't as many as i processes
        try:
            with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
                    '../../../output/batch_prune_' + str(i) + '.json')), 'r') as json_file:
                batch = json.load(json_file)
                output.update(batch)
        except Exception as e:
            print(e)
            continue
    bed_file_name = c.bed_file.split('/')[-1]
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/boundaries_prune_' + bed_file_name + '_' + str(c.ksize) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))

def run_batch(tracks, index):
    c = config.Configuration()
    for track in tracks:
        print(colorama.Fore.GREEN + '========================================================')
        print(colorama.Fore.GREEN + 'track: ', track, '@', index)
        print(colorama.Fore.GREEN, len(track))
        tracks[track] = prune_boundary_candidates(tracks[track], index)
        print(colorama.Fore.GREEN, len(track))
    print(colorama.Fore.GREEN, 'process ', index, ' done')
    # output manually, io redirection could get entangled with multiple client/servers
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_prune_' + str(index) + '.json')), 'w') as json_file:
        json.dump(tracks, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    exit()

def prune_boundary_candidates(track, index):
    # remove those candidates with high number of kmers ocurring in reference
    remove = {}
    for candidate in track:
        # skip the json key holding the number of candidates
        if candidate.find('candidates') != -1:
            continue
        kmers = track[candidate]['kmers']
        # quickly dismiess
        if not has_novel_kmers(kmers, index):
            remove[candidate] = True
            continue
        # do the more time-consuming check
        if not has_unique_novel_kmers(track, candidate, kmers, index):
            remove[candidate] = True
            continue
    for candidate in remove:
        track.pop(candidate, None)
    track['candidates'] = len(track) - 1
    return track

def has_novel_kmers(kmers, index):
    # checks if this candidate has a kmer that has not occured in the reference genome
    for kmer in kmers:
        count = count_server.get_kmer_count(kmer, index, True)
        if count == 0:
            # it is a novel kmer
            return True
    return False

def has_unique_novel_kmers(track, candidate, kmers, index):
    # checks if this candidate has a kmer that hasn't occurred in any other candidate
    for kmer in kmers:
        for break_point in track:
            # skip the wicked candidate count key
            if break_point.find('candidates') != -1:
                continue
            if candidate != break_point:
                if kmer in track[break_point]['kmers']:
                    return False
    return True

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference")
    args = parser.parse_args()
    # 
    config.configure(reference_genome = args.reference)
    #
    execute()
