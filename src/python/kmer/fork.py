import io
import os
import pwd
import sys
import json
import time

from . import (
    bed,
    sets,
    config,
    reference,
    counttable,
)

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# Classification
# ============================================================================================================================ #

def classify_sample():
    c = config.Configuration()
    #
    sample_counttable = counttable.import_sample_counttable()
    tracks = bed.import_tracks()
    #
    print(colorama.Fore.WHITE + 'classifying ... ')
    for track in tracks:
        print('--------------------------------------------------------')
        print('track: ', str(track).strip())
        track = tracks[track]
        # 
        reference = track['reference']
        variation = track['variation']
        #
        print('reference segment: ', len(reference['kmers']), ' @ ',\
            reference['head'], ' ... ', reference['tail'])
        print('variation segment: ', len(variation['kmers']), ' @ ',\
            variation['head'], ' ... ', variation['tail'])
        # 
        intersection = sets.calc_dictionary_intersection(reference['kmers'], variation['kmers'])
        print('reference/variation intersection: ', len(intersection), ' kmers')
        if intersection:
            print(intersection)
        #
        reference_score = len(calc_similarity_score(
            sets.calc_dictionary_difference(reference['kmers'], variation['kmers']), sample_counttable))
        print('fastq/reference similarity: ', reference_score)
        variation_score = len(calc_similarity_score(
            sets.calc_dictionary_difference(variation['kmers'], reference['kmers']), sample_counttable))
        print('fastq/variation similarity: ', variation_score)
        print('========================================')
        print('decision: ', ('reference' if reference_score > variation_score else
                    ('variation' if reference_score < variation_score else 'undecisive')))

def calc_similarity_score(kmers, counttable):
    print('========================================')
    result = {}
    for kmer in kmers:
        if counttable.get_kmer_counts(kmer)[0] != 0 :
            print(kmer, 'sample: ', '{:04d}'.format(counttable.get_kmer_counts(kmer)[0]))
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
        pid = os.fork()
        if pid == 0 :
            # child
            bed.export_tracks()
        else:
            print(colorama.Fore.WHITE + 'waiting for variation boundary kmers ... ')
            os.waitpid(pid, 0)
            print(colorama.Fore.WHITE + 'kmers exported.')
            classify_sample()

def main():
    config.configure()
    execute()

if __name__ == '__main__':
    main()

# ============================================================================================================================ #
# And they lived together happily forever after ...
# ============================================================================================================================ #