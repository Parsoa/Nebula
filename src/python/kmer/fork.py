import io
import os
import pwd
import sys
import json
import time

from . import (
    sets,
    config,
    reference
)

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# Constants
# ============================================================================================================================ #

class bcolors:
    COUNTTABLE = '\033[95m'
    BOUNDARIES = '\033[92m'

# ============================================================================================================================ #
# Decorators
# ============================================================================================================================ #

def measure_time(f):
    def wrapper():
        start = time.clock()
        result = f()
        end = time.clock()
        print(colorama.Fore.YELLOW + 'took ', end - start)
        return result
    return wrapper

# ============================================================================================================================ #
# CountTable
# ============================================================================================================================ #

@measure_time
def export_sample_counttable():
    c = config.Configuration()
    # 
    print(colorama.Fore.BLUE + 'searching for cached counttable ...')
    cache = c.fastq_file + '.ct'
    if os.path.isfile(cache):
        print(colorama.Fore.BLUE + 'found at ', cache)
        return
    #
    print(colorama.Fore.BLUE + 'not found, generating counttable ' + cache)
    counttable, nkmers = count_kmers_from_file(c.fastq_file)
    counttable.save(cache)
    print(colorama.Fore.BLUE + 'done')
    #
    print(colorama.Fore.BLUE + 'sample counttable cached\n', 'kmers: ', nkmers,\
        '\nsize: ', os.stat(cache).st_size)
    return

@measure_time
def import_sample_counttable():
    print(colorama.Fore.MAGENTA + 'importing counttable ...')
    c = config.Configuration()
    cache = c.fastq_file + '.ct'
    counttable = khmer.Counttable.load(cache)
    print(colorama.Fore.MAGENTA + 'done')
    return counttable

def count_kmers_from_file(file):
    c = config.Configuration()
    #
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counttable.consume_seqfile(file)
    #
    return counttable, nkmers

# ============================================================================================================================ #
# SV Tracks
# ============================================================================================================================ #

def export_tracks():
    c = config.Configuration()
    bedtools = pybedtools.BedTool(c.bed_file)
    #
    kmers = {}
    for track in bedtools:
        print(colorama.Fore.GREEN + '--------------------------------------------------------')
        #
        key = track.chrom.strip() + '_' + str(track.start).strip()  + '_'+ str(track.end).strip()
        print(colorama.Fore.GREEN + 'track: ', key)
        # 
        reference_boundaries, variation_boundaries = extract_track_boundaries(track)
        reference_kmers = count_boundary_kmers(reference_boundaries)
        variation_kmers = count_boundary_kmers(variation_boundaries)
        #
        kmers[key] = {
            'reference':{
                'head':  reference_boundaries['head'],
                'tail':  reference_boundaries['tail'],
                'kmers': reference_kmers
            },
            'variation':{
                'head':  variation_boundaries['head'],
                'tail':  variation_boundaries['tail'],
                'kmers': variation_kmers
            }
        }
    # print(colorama.Fore.GREEN, kmers)
    path = os.path.join(c.output_directory, 'boundaries.json')
    print(colorama.Fore.GREEN + 'exporting tracks to ', path)
    with open(path, 'w') as json_file:
        json.dump(kmers, json_file, sort_keys=True, indent=4, separators=(',', ': '))

def import_tracks():
    c = config.Configuration()
    with io.open(os.path.join(c.output_directory, 'boundaries.json'), 'r') as json_file:
        tracks = json.load(json_file)
        return tracks

def extract_track_boundaries(track):
    c = config.Configuration()
    interval = pybedtools.Interval(chrom = track.chrom, start = track.start - c.ksize,\
        end = track.end + c.ksize)
    bedtool = pybedtools.BedTool(str(interval), from_string = True)
    f = open((bedtool.sequence(fi = reference.ReferenceGenome().fasta)).seqfn)
    for i, line in enumerate(f):
        line = line.strip()
        if i == 1:
            # print('#', line, '#')
            head = line[0:2 * c.ksize]
            tail = line[-2 * c.ksize:]
            # print('head: ', head)
            # print('tail: ', tail)
            # inverse this sequence
            line = line[::-1]
            # print('#', line, '#')
            inverse_head = line[0:2 * c.ksize]
            inverse_tail = line[-2 * c.ksize:]
            return {'head': head.upper(), 'tail': tail.upper()}, {'head': inverse_head.upper(), 'tail': inverse_tail.upper()}

def count_boundary_kmers(boundaries):
    c = config.Configuration()
    #
    kmers = {}
    kmers = count_kmers(boundaries['head'], c.ksize, kmers)
    kmers = count_kmers(boundaries['tail'], c.ksize, kmers)
    print(kmers)
    #
    return kmers

def count_kmers(str, k, kmers):
    for i in range(0, len(str) - k):
        kmer = str[i:i + k]
        # print(kmer, ' ',  len(kmer))
        if not kmer in kmers :
            kmers[kmer] = 1
        else :
            kmers[kmer] = kmers[kmer] + 1
    return kmers

# ============================================================================================================================ #
# Reference CountTable
# ============================================================================================================================ #

# TODO

# ============================================================================================================================ #
# Classification
# ============================================================================================================================ #

def classify_sample():
    c = config.Configuration()
    #
    sample_counttable = import_sample_counttable()
    tracks = import_tracks()
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
# Configuration
# ============================================================================================================================ #

def configure():
    if sys.platform == "darwin":
        print('Running on Mac OS X')
        khmer_table_size = 16e7
        khmer_num_tables = 4
        reference.ReferenceGenome(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/hg38.fa')))
        fastq_file = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1.samtoolsversion.head.tiny.fq'))
    else:
        print('Running on Linux')
        khmer_table_size = 5e9
        khmer_num_tables = 4
        fastq_file = '/share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq'
        reference.ReferenceGenome('/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa')
    config.Configuration(
        ksize = 25,
        khmer_table_size = khmer_table_size,
        khmer_num_tables = khmer_num_tables,
        fastq_file = fastq_file,
        bed_file = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1.inversions_hg38.bed')),
        output_directory = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../output'))
    )
    colorama.init()

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

if __name__ == '__main__':
    configure()
    pid = os.fork()
    if pid == 0:
        # child
        export_sample_counttable()
    else:
        print(colorama.Fore.WHITE + 'waiting for couttable creation ... ')
        os.waitpid(pid, 0)
        print(colorama.Fore.WHITE + 'counttable exported.')
        pid = os.fork()
        if pid == 0 :
            # child
            export_tracks()
        else:
            print(colorama.Fore.WHITE + 'waiting for variation boundary kmers ... ')
            os.waitpid(pid, 0)
            print(colorama.Fore.WHITE + 'kmers exported.')
            classify_sample()

# ============================================================================================================================ #
# And they lived together happily forever after ...
# ============================================================================================================================ #