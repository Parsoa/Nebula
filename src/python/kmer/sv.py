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
    sample_counttable = counttable.import_sample_counttable()
    #
    for track in bedtools:
        find_track_boundaries(track, sample_counttable)

@commons.measure_time
def find_track_boundaries(track, counttable):
    print(colorama.Fore.GREEN + '========================================================')
    print(colorama.Fore.GREEN + 'track: ', str(track).strip())
    #
    max = None
    for start in range(-radius, radius) :
        for end in range(-radius, radius) :
            try :
                boundary = (start, end)
                score = calc_boundary_score(track, boundary, counttable)
                print(colorama.Fore.GREEN + 'score: ', score)
                max = (score, boundary) if not max else\
                    (score, boundary) if score > max[0] else max
            except Exception as e:
                print(e)
                print(colorama.Fore.RED + 'boundary error, skipping')
    print(colorama.Fore.GREEN + '########################################################')
    print('choice: ', max)

# @commons.measure_time
def calc_boundary_score(t, boundary, counttable):
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
    score = len(calc_similarity_score(variation_kmers, counttable))
    return score

def calc_similarity_score(kmers, counttable):
    result = {}
    for kmer in kmers:
        if counttable.get_kmer_counts(kmer)[0] != 0 :
            print(kmer, '{:04d}'.format(counttable.get_kmer_counts(kmer)[0]))
            result[kmer] = True
    return result

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

def execute():
    pid = os.fork()
    if pid == 0:
        # child
        counttable.export_sample_counttable()
    else:
        print(colorama.Fore.WHITE + 'waiting for couttable creation ... ')
        os.waitpid(pid, 0)
        print(colorama.Fore.WHITE + 'counttable exported.')
        refine_variation_boundaries()

if __name__ == '__main__':
    config.configure()
    execute()
