import io
import os
import re
import pwd
import sys
import copy
import json
import time
import atexit
import argparse
import traceback

from kmer import (
    config,
    commons,
)

import colorama

# ============================================================================================================================ #
# Global variables, no object hierarchy here everything has to be here
# ============================================================================================================================ #

kmers = {}
job_name = "CountKmersExactJob_"
num_threads = 12
thread_pool = None
previous_job_name = "NovelKmersJob_"

# ============================================================================================================================ #
# filesystem helpers
# ============================================================================================================================ #

def get_output_directory():
    c = config.Configuration()
    bed_file_name = c.bed_file.split('/')[-1]
    return os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../output/' + bed_file_name + '/' + str(c.ksize) + '/'))

def get_previous_job_directory():
    # get rid of the final _
    return os.path.abspath(os.path.join(get_output_directory(), previous_job_name[:-1]))

def get_current_job_directory():
    # get rid of the final _
    return os.path.abspath(os.path.join(get_output_directory(), job_name[:-1]))

def create_output_directories():
    dir = get_output_directory()
    if not os.path.exists(dir):
        os.makedirs(dir)
    dir = get_current_job_directory()
    if not os.path.exists(dir):
        os.makedirs(dir)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

def parse_fastq(index, fastq_file):
    name = None
    HEADER_LINE = 0
    SEQUENCE_LINE = 1
    THIRD_LINE = 2
    QUALITY_LINE = 3
    state = HEADER_LINE
    # need to skip invalid lines
    line = fastq_file.readline()
    ahead = fastq_file.readline()
    n = 0
    m = 0
    t = time.time()
    # TODO: to make IO waits more meaningful, load a 1000 lines or so per iteration
    while ahead:
        if state == HEADER_LINE:
            if line[0] == '@' and ahead[0] != '@':
                if fastq_file.tell() >= (index + 1) * fastq_file_chunk_size:
                    print(index, 'reached segment boundary')
                    break
                state = SEQUENCE_LINE
                name = line[:-1] # ignore the EOL character
        elif state == SEQUENCE_LINE:
            state = THIRD_LINE
            seq = line[:-1] # ignore the EOL character
            n += 1
            if n == 100000:
                n = 0
                m += 1
                c = fastq_file.tell() - index * fastq_file_chunk_size
                s = time.time()
                p = c / float(fastq_file_chunk_size)
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                print(index, 'progress:', p, 'took: ', s - t, 'ETA: ', e)
            yield seq, name
        elif state == THIRD_LINE:
            state = QUALITY_LINE
        elif state == QUALITY_LINE:
            state = HEADER_LINE
        line = ahead
        ahead = fastq_file.readline()
    print(index, ' end of input')

def get_all_kmers(read, k):
    kmers = []
    for i in range(0, len(read) - k + 1):
        kmer = read[i : i + k]
        kmers.append(kmer)
    return kmers

def find_thread_count():
    c = config.Configuration()
    max_index = 0
    for index in range(0, c.max_threads):
        path = os.path.join(get_previous_job_directory(), 'batch_' + str(index) + '.json')
        if os.path.isfile(path):
            max_index = index + 1
    global num_threads
    num_threads = max_index
    global thread_pool
    thread_pool = ThreadPool(num_threads) 

def load_inputs():
    c = config.Configuration()
    #
    global fastq_file_chunk_size
    fastq_file_chunk_size = math.ceil(os.path.getsize(fastq_file.name) / float(num_threads))
    # 
    path = os.path.join(get_previous_job_directory(), 'merge.json')
    tracks = {}
    with open(path, 'r') as json_file:
        batch = json.load(json_file)
        for track in batch:
            tracks[track] = batch[track]
    #
    global kmers
    for track in tracks:
        novel_kmers = tracks[track]['novel_kmers']
        for kmer in novel_kmers:
            kmers[kmer] = 0

def transform(index):
    c = config.Configuration()
    fastq_file = open(c.fastq_file, 'r')
    fastq_file.seek(index * fastq_file_chunk_size, 0)
    for read, name in parse_fastq(index, fastq_file):
        for kmer in get_all_kmers(read, c.ksize):
            canon = get_canonical_kmer_representation(kmer)
            global kmers
            if canon in kmers: 
                kmers[canon] += 1

def join():
    c = config.Configuration()
    # no need to merge kmer counts
    # reindex based on track
    path = os.path.join(.get_previous_job_directory(), 'merge.json')
    tracks = {}
    with open(path, 'r') as json_file:
        batch = json.load(json_file)
        for track in batch:
            tracks[track] = batch[track]
    for track in tracks:
        novel_kmers = tracks[track]['novel_kmers']
        for novel_kmer in novel_kmers:
            novel_kmers[novel_kmer]['actual_count'] = kmers[novel_kmer]
    # 
    with open(os.path.join(get_current_job_directory(), 'merge.json'), 'w') as json_file:
        json.dump(tracks, json_file, sort_keys = True, indent = 4, separators = (',', ': '))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

def execute():
    find_thread_count()
    create_output_directories()
    global kmers
    kmers = load_inputs()
    thread_pool.map(transform, range(0, num_threads))
    thread_pool.close()
    thread_pool.join()
    join()
