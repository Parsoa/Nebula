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
import pybedtools

# ============================================================================================================================ #
# Decorators
# ============================================================================================================================ #

def measure_time(f):
    def wrapper():
        start = time.clock()
        result = f()
        end = time.clock()
        print('took ', end - start)
        return result
    return wrapper

def conditionally(f):
    print(__name__)
    if __name__ == '__main__':
        print('not profiling')
        return f
    print('profiling')
    return profile(f)

# ============================================================================================================================ #
# CountTable
# ============================================================================================================================ #

@measure_time
def export_sample_counttable():
    c = config.Configuration()
    # 
    print('searching for cached counttable ...')
    cache = c.fastq_file + '.ct'
    if os.path.isfile(cache):
        print('found at ', cache)
        print('done')
        return
    #
    print('not found, generating counttable ' + cache)
    counttable, nkmers = count_kmers_from_file(c.fastq_file)
    counttable.save(cache)
    print('done')
    #
    print('sample counttable cached\n', 'kmers: ', nkmers,\
        '\nsize: ', os.stat(cache).st_size)
    return

@measure_time
def import_sample_counttable():
    print('importing counttable ...')
    c = config.Configuration()
    cache = c.fastq_file + '.ct'
    return khmer.Counttable.load(cache)

# ============================================================================================================================ #
# CountTable
# ============================================================================================================================ #

def export_tracks():
    c = config.Configuration()
    bedtools = pybedtools.BedTool(c.bed_file)
    #
    kmers = {}
    for track in bedtools:
        print('--------------------------------------------------------')
        print('track: ', str(track).strip())
        #
        key = track.chrom.strip() + '_' + str(track.start).strip()  + '_'+ str(track.end).strip()
        # 
        reference_boundaries, variation_boundaries = extract_track_boundaries(track)
        reference_kmers = count_boundary_kmers(reference_boundaries)
        variation_kmers = count_boundary_kmers(variation_boundaries)
        #
        kmers[key] = {'reference': reference_kmers, 'variation': variation_kmers}
    json_file = open(os.path.join(c.output_directory, 'boundaries.json'), 'w+')
    json.dump(obj, json_file)

def import_tracks():
    c = config.Configuration()
    with io.open(os.path.join(c.output_directory, 'boundaries.json'), 'w') as json_file:
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
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    #
    kmers = {}
    #
    counttable.consume(boundaries['head'])
    counttable.consume(boundaries['tail'])
    #
    for sequence in [boundaries['head'], boundaries['tail']]:
        for kmer in counttable.get_kmers(sequence) :
            kmers[kmer] = counttable.get_kmer_counts(kmer)
    return kmers

# ============================================================================================================================ #
# Reference CountTable
# ============================================================================================================================ #

@measure_time
def count_reference_kmers():
    print('reading reference genome ... ')
    c = config.Configuration()
    #
    counttable, nkmers = count_kmers_from_file(reference.ReferenceGenome().path)
    #
    print('reference counttable cached: ', nkmers, ' kmers recorded')
    return counttable

def count_kmers_from_file(file):
    c = config.Configuration()
    #
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counttable.consume_seqfile(file)
    #
    return counttable, nkmers

# ============================================================================================================================ #
# Classification
# ============================================================================================================================ #

@conditionally
def classify_sample():
    c = config.Configuration()
    #
    sample_counttable = import_sample_counttable()
    tracks = import_tracks()
    #
    for track in tracks:
        print('--------------------------------------------------------')
        print('track: ', str(track).strip())
        # 
        reference_kmers = track['reference']
        variation_kmers = track['variation']
        #
        print('reference segment: ', len(reference_kmers.keys()), ' @ ',\
            reference_boundaries['head'], ' ... ', reference_boundaries['tail'])
        print('variation segment: ', len(variation_kmers.keys()), ' @ ',\
            variation_boundaries['head'], ' ... ', variation_boundaries['tail'])
        # 
        intersection = sets.calc_dictionary_intersection(reference_kmers, variation_kmers)
        print('reference/variation intersection: ', len(intersection), ' kmers')
        if intersection:
            print(intersection)
        #
        reference_score = len(calc_similarity_score(
            sets.calc_dictionary_difference(reference_kmers, variation_kmers), sample_counttable))
        print('fastq/reference similarity: ', reference_score)
        variation_score = len(calc_similarity_score(
            sets.calc_dictionary_difference(variation_kmers, reference_kmers), sample_counttable))
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
    # print('cwd: {}'.format(os.getcwd()))
    khmer_table_size = 16e9
    khmer_num_tables = 4
    if sys.platform == "darwin":
        print('Running on Mac OS X')
        reference.ReferenceGenome(os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/hg38.fa')))
        fastq_file = os.path.abspath(os.path.join(os.path.dirname(__file__),\
            '../../../data/CHM1.samtoolsversion.head.tiny.fq'))
    else:
        print('Running on Linux')
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

# ============================================================================================================================ #
# Execution
# ============================================================================================================================ #

if __name__ == '__main__':
    configure()
    pid = os.fork()
    if pid == 0:
        print('child process')
        # export_sample_counttable()
    else:
        print('root process, waiting for couttable creation ... ')
        os.waitpid(pid, 0)
        print('root process, counttable exported.')
        print('root process, extracting reference/variation kmers')
        pid = os.fork()
        if pid == 0 :
            print('child process, extracting track boundaries')
            export_tracks()
            print('child process, done')
        else:
            print('root process, waiting for kmers ... ')
            os.waitpid(pid, 0)
            print('root process, kmers exported.')
            print('proceeding ... ')
            classify_sample()

 