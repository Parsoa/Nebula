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
import plotly.offline as plotly
import plotly.graph_objs as go

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

def execute():
    c = config.Configuration()
    # find the number of processes to spawn
    max_index = 0
    for index in range(0, c.num_threads):
        path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_' + str(index) + '.json'))
        if os.path.isfile(path):
            max_index = index + 1
    # run each batch
    children = {}
    for index in range(0, max_index):
        pid = os.fork()
        if pid == 0:
            # forked process
            path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_' + str(index) + '.json'))
            with open(path, 'r') as json_file:
                print('reading batch ', index)
                batch = json.load(json_file)
                run_batch(batch, index)
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
    draw_distribution_charts(output)
    clean_up()

def clean_up():
    for i in range(0, c.num_threads):
        # might fail because there weren't as many as i processes
        path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_prune_' + str(i) + '.json'))
        if os.path.isfile(path):
                os.remove(path)

def draw_distribution_charts(tracks):
    c = config.Configuration()
    bins = {}
    for track in tracks:
        n = tracks[track]['breakpoint_without_novel']
        if not n in bins:
            bins[n] = 1
        else:
            bins[n] = bins[n] + 1
    data = [go.Bar(
        x = bins.keys(),
        y = list(map(lambda x: bins[x], bins.keys()))
    )]
    bed_file_name = c.bed_file.split('/')[-1]
    plotly.plot(data, filename = os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/boundaries_prune_' + bed_file_name + '_' + str(c.ksize) + '.html')))

def run_batch(tracks, index):
    c = config.Configuration()
    for track in tracks:
        print(colorama.Fore.GREEN + '========================================================')
        print(colorama.Fore.GREEN + 'track: ', track, '@', index)
        tracks[track] = aggregate_novel_kmers(tracks[track], index)
    print(colorama.Fore.GREEN, 'process ', index, ' done')
    # output manually, io redirection could get entangled with multiple client/servers
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_prune_' + str(index) + '.json')), 'w') as json_file:
        json.dump(tracks, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    exit()

def aggregate_novel_kmers(track, index):
    n = 0
    remove = {}
    for candidate in track:
        # skip the json key holding the number of candidates
        if candidate.find('candidates') != -1:
            continue
        kmers = track[candidate]['kmers']
        remove[candidate] = True
        if not has_novel_kmers(kmers, index):
            n += 1
        # novel_kmers = get_novel_kmers(kmers, index)
        # track[candidate]['novel_kmers'] = novel_kmers
        # track[candidate]['novel_kmer_count'] = len(novel_kmers)
        # novel_kmer_count += len(novel_kmers)
        # n = n + 1
    # remove those break points without uniquely novel kmers
    for candidate in remove:
        track.pop(candidate, None)
    # cleanup unwanted keys
    # for candidate in track:
    #     # skip the json key holding the number of candidates
    #     if candidate.find('candidates') != -1:
    #         continue
    #     track[candidate].pop('kmers', None)
    #     track[candidate].pop('reference_kmers', None)
    track['breakpoint_without_novel'] = n
    # track['average_novel_kmer_count'] = -1 if n < 1 else novel_kmer_count / n
    # track['candidates'] = n
    return track

def get_novel_kmers(kmers, index):
    novel_kmers = {}
    for kmer in kmers:
        if not kmer in novel_kmers:
            count = count_server.get_kmer_count(kmer, index, True)
            if count != 0:
                # this is a novel kmer
                novel_kmers[kmer] = True
    return novel_kmers

def prune_boundary_candidates(track, index):
    # remove those candidates with high number of kmers ocurring in reference
    remove = {}
    contigs = {}
    for candidate in track:
        # skip the json key holding the number of candidates
        if candidate.find('candidates') != -1:
            continue
        kmers = track[candidate]['kmers']
        contig = track[candidate]['boundary']
        if not contig in contigs:
            contigs[contig] = 1
        else:
            contigs[contig] = contigs[contig] + 1
        # quickly dismiess
        if not has_novel_kmers(kmers, index):
           remove[candidate] = True
           continue
        # do the more time-consuming check
        if not has_unique_novel_kmers(track, candidate, kmers, index):
           remove[candidate] = True
           continue
        track[candidate].pop('reference_kmers', None)
    for candidate in remove:
        track.pop(candidate, None)
    track['contig_count'] = len(contigs)
    track['candidates'] = len(track) - 1
    track['contigs'] = contigs
    return track

def has_novel_kmers(kmers, index):
    # checks if this candidate has a kmer that has not occured in the reference genome
    for kmer in kmers:
        if is_kmer_novel(kmer, index):
            return True
    return False

def is_kmer_novel(kmer, index):
    count = count_server.get_kmer_count(kmer, index, True)
    return count == 0

def has_unique_novel_kmers(track, candidate, kmers, index):
    # checks if this candidate has a novel kmer that hasn't occurred in any other candidate
    for kmer in kmers:
        if is_kmer_novel(kmer, index):
            found = False
            for break_point in track:
                # skip the wicked candidate count key
                if break_point.find('candidates') != -1:
                    continue
                if candidate != break_point:
                    # this kmer appears in at least one other break point so no need to go further
                    if kmer in track[break_point]['kmers']:
                        found = True
                        break
            # we didn't find this novel kmer anywhere so it should be unique, no need to check others
            if not found:
                return True
    # we haven't returned yet so we didn't find any uniquely novel kmers
    return False

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