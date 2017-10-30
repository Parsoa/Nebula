import io
import os
import pwd
import sys
import json
import time

from . import (
    sets,
    config,
    commons,
    reference,
)

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# BED Tracks
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
            'start': track.start,
            'end'  : track.end,
            'reference': {
                'head':  reference_boundaries['head'],
                'tail':  reference_boundaries['tail'],
                'kmers': reference_kmers
            },
            'variation': {
                'head':  variation_boundaries['head'],
                'tail':  variation_boundaries['tail'],
                'kmers': variation_kmers
            }
        }
    # print(colorama.Fore.GREEN, kmers)
    path = os.path.join(c.output_directory, 'boundaries.json')
    print(colorama.Fore.GREEN + 'exporting tracks to ', path)
    with open(path, 'w') as json_file:
        json.dump(kmers, json_file, indent = 4, separators = (',', ': '))

def import_tracks():
    c = config.Configuration()
    with io.open(os.path.join(c.output_directory, 'boundaries.json'), 'r') as json_file:
        tracks = json.load(json_file)
        return tracks

def extract_sequence(track):
    c = config.Configuration()
    interval = pybedtools.Interval(chrom = track.chrom, start = track.start,\
        end = track.end)
    bedtool = pybedtools.BedTool(str(interval), from_string = True)
    f = open((bedtool.sequence(fi = reference.ReferenceGenome().fasta)).seqfn)
    for i, line in enumerate(f):
        if i == 1:
            line = line.strip().upper()
            return line

def extract_track_boundaries(track):
    c = config.Configuration()
    interval = pybedtools.Interval(chrom = track.chrom, start = track.start - c.ksize,\
        end = track.end + c.ksize)
    bedtool = pybedtools.BedTool(str(interval), from_string = True)
    f = open((bedtool.sequence(fi = reference.ReferenceGenome().fasta)).seqfn)
    for i, line in enumerate(f):
        if i == 1:
            line = line.strip().upper()
            # print(colorama.Fore.WHITE + line[:c.ksize], colorama.Fore.BLUE + (line[c.ksize : -c.ksize]), colorama.Fore.WHITE + line[-c.ksize:])
            head = line[0:2 * c.ksize]
            tail = line[-2 * c.ksize:]
            # reverse-complement this sequence
            # print(colorama.Fore.WHITE + "sequence: ", line[:c.ksize], '#',\
                #‌ colorama.Fore.BLUE + complement_sequence((line[c.ksize : -c.ksize])[::-1]), '#',\
                # colorama.Fore.WHITE + line[-c.ksize:])
            # line = line[:c.ksize] + complement_sequence((line[c.ksize : -c.ksize])[::-1]) + line[-c.ksize:]
            inverse_head = line[0:2 * c.ksize]
            inverse_tail = line[-2 * c.ksize:]
            # print(colorama.Fore.WHITE + "boundary: ", colorama.Fore.BLUE + inverse_head, '...', inverse_tail)
            return {'head': head, 'tail': tail}, {'head': inverse_head, 'tail': inverse_tail}

def complement_sequence(seq):
    # A-> C and C->A
    seq = seq.replace('A', 'Z')
    seq = seq.replace('T', 'A')
    seq = seq.replace('Z', 'T')
    #
    seq = seq.replace('G', 'Z')
    seq = seq.replace('C', 'G')
    seq = seq.replace('Z', 'C')
    #
    return seq