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
import plotly.graph_objs as graph_objs

events = {}

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

def get_novel_kmers(kmers, index):
    novel_kmers = {}
    for kmer in kmers:
        if not kmer in novel_kmers:
            count = count_server.get_kmer_count(kmer, index, True)
            if count != 0:
                # this is a novel kmer
                novel_kmers[kmer] = True
    return novel_kmers

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
# map functions
# ============================================================================================================================ #

def prune_boundary_candidates(track, track_name, index):
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

# ============================================================================================================================ #

def map_novel_kmer_overlap(track, track_name, index):
    novel_kmers = track['novel_kmers']
    overlap = {}
    global events
    for event in events:
        if event != track_name:
            for kmer in novel_kmers:
                if kmer in events[event]['novel_kmers']:
                    if not kmer in overlap:
                        overlap[kmer] = []
                    overlap[kmer].append(event)
    track['overlap'] = overlap
    track['overlap_percentage'] = -1 if len(novel_kmers) == 0 else float(len(overlap)) / float(len(novel_kmers))
    track.pop('novel_kmers', None)
    return track

def draw_novel_kmer_overlap_plot(tracks, job_name):
    path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/merge_' + job_name + '.html'))
    x = list(map(lambda x: tracks[x]['overlap_percentage'], tracks))
    trace = graph_objs.Histogram(
        x = x,
        histnorm = 'count',
        xbins = dict(
            start = 0.0,
            end = 1.0,
            size = 0.05
        )
    )
    layout = graph_objs.Layout(
        title = 'Novel kmer Overlap'
    )
    fig = graph_objs.Figure(data = [trace], layout = layout)
    plotly.plot(fig, filename = path)

# ============================================================================================================================ #

def aggregate_novel_kmers(track, track_name, index):
    remove = {}
    novel_kmers = {}
    for candidate in track:
        # remove all candidates eventually
        remove[candidate] = True
        # skip the json key holding the number of candidates
        if candidate.find('candidates') != -1:
            continue
        kmers = track[candidate]['kmers']
        for kmer in kmers:
            if is_kmer_novel(kmer, index):
                if not kmer in novel_kmers:
                    novel_kmers[kmer] = 1
                else:
                    novel_kmers[kmer] = novel_kmers[kmer] + 1
    # remove every candidate, this has to be lightweight
    for candidate in remove:
        track.pop(candidate, None)
    # keep only those with high enough coverage
    novel_kmers = list(filter(lambda kmer: novel_kmers[kmer] > 5, novel_kmers))
    track['novel_kmers'] = novel_kmers
    return track

# ============================================================================================================================ #
# concurrency helpers
# ============================================================================================================================ #

def find_thread_count():
    c = config.Configuration()
    max_index = 0
    for index in range(0, c.max_threads):
        path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_' + str(index) + '.json'))
        if os.path.isfile(path):
            max_index = index + 1
    return max_index

def distribute_workload(job_name, previous_job_name, func, num_threads):
    children = {}
    for index in range(0, num_threads):
        pid = os.fork()
        if pid == 0:
            # forked process
            path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_' + previous_job_name + str(index) + '.json'))
            with open(path, 'r') as json_file:
                print('reading batch ', index)
                batch = json.load(json_file)
                run_batch(batch, index, func, job_name)
        else:
            # main process
            children[pid] = True
            print('spawned child ', pid)
    return children

def run_batch(batch, index, func, job_name):
    c = config.Configuration()
    for track in batch:
        # print(colorama.Fore.GREEN + '========================================================')
        # print(colorama.Fore.GREEN + 'batch: ', batch, '@', index)
        batch[track] = func(batch[track], track, index)
    print(colorama.Fore.GREEN, 'process ', index, ' done')
    # output manually, io redirection could get entangled with multiple client/servers
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_' + job_name + str(index) + '.json')), 'w') as json_file:
        json.dump(batch, json_file, sort_keys=True, indent=4, separators=(',', ': '))
    exit()

def wait_for_children(children):
    while True:
        (pid, e) = os.wait()
        children.pop(pid, None)
        print(colorama.Fore.RED, 'pid ', pid, 'finished')
        if len(children) == 0:
            break

def merge_outputs(job_name, num_threads, merge_function):
    c = config.Configuration()
    output = {}
    for i in range(0, num_threads):
        with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
                '../../../output/batch_' + job_name + str(i) + '.json')), 'r') as json_file:
            batch = json.load(json_file)
            output.update(batch)
    if merge_function:
        merge_function(output, job_name)
    bed_file_name = c.bed_file.split('/')[-1]
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/merge_' + job_name + bed_file_name + '_' + str(c.ksize) + '.json')), 'w') as json_file:
        json.dump(output, json_file, sort_keys=True, indent=4, separators=(',', ': '))

def clean_up(job_name, num_threads):
    c = config.Configuration()
    for i in range(0, num_threads):
        # might fail because there weren't as many as i processes
        path = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output/batch_' + job_name + str(i) + '.json'))
        os.remove(path)

# ============================================================================================================================ #
# job descriptions
# ============================================================================================================================ #

def find_high_coverage_novel_kmers():
    c = config.Configuration()
    previous_job_name = ''
    job_name = 'novel_'
    # 
    num_threads = find_thread_count()
    children = distribute_workload(job_name, previous_job_name, aggregate_novel_kmers, num_threads)
    wait_for_children(children)
    print('all forks done, merging output ...')
    merge_outputs(job_name, num_threads, None)
    # clean_up(job_name, num_threads)

def find_novel_kmer_overlap_map():
    c = config.Configuration()
    previous_job_name = 'novel_'
    job_name = 'overlap_'
    # prepare
    bed_file_name = c.bed_file.split('/')[-1]
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/merge_' + previous_job_name + bed_file_name + '_' + str(c.ksize) + '.json')), 'r') as json_file:
        global events
        events = json.load(json_file)
    # 
    num_threads = find_thread_count()
    children = distribute_workload(job_name, previous_job_name, map_novel_kmer_overlap, num_threads)
    wait_for_children(children)
    print('all forks done, merging output ...')
    merge_outputs(job_name, num_threads, draw_novel_kmer_overlap_plot)
    # clean_up(job_name, num_threads)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #            

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
    find_novel_kmer_overlap_map()
    # find_high_coverage_novel_kmers()