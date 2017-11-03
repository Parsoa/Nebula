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
        with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_prune_' + str(i) + '.json')), 'r') as json_file:
            batch = json.load(json_file)
            output.update(batch)
    bed_file_name = c.bed_file.split('/')[-1]
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/boundaries_' + bed_file_name + '_' + str(c.ksize) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))

def run_batch(tracks, index):
    c = config.Configuration()
    output = {}
    for track in tracks:
        print(colorama.Fore.GREEN + '========================================================')
        print(colorama.Fore.GREEN + 'track: ', track, '@', index)
        tracks[track] = prune_boundary_candidates(tracks[track], index)
    print(colorama.Fore.GREEN, 'process ', index, ' done')
    # output manually, io redirection could get entangled with multiple client/servers
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_prune_' + str(index) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    exit()

def prune_boundary_candidates(track, index):
    print(track)
    # remove those candidates with high number of kmers ocurring in reference
    remove = {}
    for candidate in track:
        kmers = candidate['kmers']
        prune_kmers(kmers)
        if len(kmers) == 0:
            remove[candidate] = True
            continue
        candidate['kmers'] = kmers
    for candidate in remove:
        tracks.pop(candidate, None)
    return track

def prune_kmers(kmers):
    remove = {}
    for kmer in kmers:
        count = count_server.get_kmer_count(kmer, index, True)
        if count:
            remove.append(kmer)
    for kmer in remove:
        kmers.pop(kmer, None)

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference")
    args = parser.parse_args()
    # 
    config.configure(args)
    # count_server.run_server('hg19')
    execute()
    count_server.kill()
