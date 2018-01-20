import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import atexit
import argparse
import traceback

from multiprocessing.dummy import Pool as ThreadPool

from kmer import (
    config,
    commons,
)

import colorama

# ============================================================================================================================ #
# Global variables, no object hierarchy here everything has to be here
# ============================================================================================================================ #

kmers = {}
job_name = "CountKmersExactThreadedJob_"
num_threads = 12
thread_pool = None
previous_job_name = "NovelKmerJob_"

# ============================================================================================================================ #
# filesystem helpers
# ============================================================================================================================ #

def get_output_directory():
    c = config.Configuration()
    bed_file_name = c.bed_file.split('/')[-1]
    return os.path.abspath(os.path.join(os.path.dirname(__file__),\
        '../../../../output/' + bed_file_name + '/' + str(c.ksize) + '/'))

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
    n = 0
    m = 0
    lines = []
    t = time.time()
    chunk_size = 1000
    for i in range(0, chunk_size):
        lines.append(fastq_file.readline())
    line = lines.pop(0)
    ahead = lines.pop(0)
    while ahead:
        if state == HEADER_LINE:
            # both the header and phred line can start with @ so we need to distinguish between them
            # this can happen the first time we are looking for a header as our offset may start from
            # a phred line
            if line[0] == '@' and ahead[0] != '@':
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
        # we have reached the end of this chunk of lines, read another chunk
        if len(lines) == 0:
            j = 0
            while True:
                if state == HEADER_LINE:
                    # we have reached a header line in the next segment of the file, stop
                    if fastq_file.tell() >= (index + 1) * fastq_file_chunk_size:
                        break
                l = fastq_file.readline()
                # end of file probably, only happens for the last segment
                if not l:
                    break
                lines.append(l)
                j = j + 1
                if j == chunk_size():
                    break
        ahead = lines.pop(0)
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

def load_inputs():
    c = config.Configuration()
    #
    global fastq_file_chunk_size
    fastq_file = open(c.fastq_file, 'r')
    fastq_file_chunk_size = math.ceil(os.path.getsize(fastq_file.name) / float(num_threads))
    #
    print('loading tracks ...') 
    path = os.path.join(get_previous_job_directory(), 'merge.json')
    tracks = {}
    with open(path, 'r') as json_file:
        batch = json.load(json_file)
        for track in batch:
            tracks[track] = batch[track]
    #
    print('re-indexing on kmers')
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
    # load the sv break points
    path = os.path.join(get_previous_job_directory(), 'merge.json')
    tracks = {}
    with open(path, 'r') as json_file:
        batch = json.load(json_file)
        for track in batch:
            tracks[track] = batch[track]
    # re-index based on track
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
    global thread_pool
    thread_pool = ThreadPool(num_threads) 
    thread_pool.map(transform, range(0, num_threads))
    thread_pool.close()
    thread_pool.join()
    join()

if __name__ == '__main__':
    config.init()
    execute()
