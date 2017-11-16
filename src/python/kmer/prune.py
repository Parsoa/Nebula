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
    map_reduce,
    break_point,
    count_server,
)
from kmer.sv import StructuralVariation, Inversion, Deletion

import khmer
import colorama
import pybedtools
import plotly.offline as plotly
import plotly.graph_objs as graph_objs

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
# ============================================================================================================================ #
# MapReduce job to extracting high-coverage novel kmers
# ============================================================================================================================ #
# ============================================================================================================================ #

class HighCoverageNovelJob(map_reduce.Job):

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
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
                if is_kmer_novel(kmer, self.index):
                    if not kmer in novel_kmers:
                        novel_kmers[kmer] = kmers[kmer]
        # remove every candidate, this has to be lightweight
        for candidate in remove:
            track.pop(candidate, None)
        # keep only those with high enough coverage
        novel_kmers = list(filter(lambda kmer: novel_kmers[kmer] > 5, novel_kmers))
        track['novel_kmers'] = novel_kmers
        return track

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job to calculate novel kmer overlap
# ============================================================================================================================ #
# ============================================================================================================================ #

class NovelKmerOverlapJob(map_reduce:Job):

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare():
        with open(os.path.join(self.get_previous_job_directory(), 'merge.json'), 'r') as json_file:
        self.events = json.load(json_file)

    def transform(self, track, track_name):
        novel_kmers = track['novel_kmers']
        overlap = {}
        for event in self.events:
            if event != track_name:
                for kmer in novel_kmers:
                    if kmer in events[event]['novel_kmers']:
                        if not kmer in overlap:
                            overlap[kmer] = []
                        overlap[kmer].append(event)
        track['overlap'] = overlap
        track['overlap_percentage'] = -1 if len(novel_kmers) == 0 else float(len(overlap)) / float(len(novel_kmers))
        return track

    def plot(self, tracks):
        self.plot_novel_kmer_overlap_count(tracks)
        self.plot_novel_kmer_overlap_percentage(tracks)

    def plot_novel_kmer_overlap_count(self, tracks):
        path = os.path.join(self.get_current_job_directory(), 'plot_novel_kmer_overlap_count.html')
        x = list(map(lambda x: len(tracks[x]['novel_kmers']) - len(track[x]['overlap']), tracks))
        trace = graph_objs.Histogram(
            x = x,
            histnorm = 'count',
            xbins = dict(
                start = 0.0,
                end = 1.0,
            )
        )
        layout = graph_objs.Layout(title = 'Non-Overlapping Novel kmer Count')
        fig = graph_objs.Figure(data = [trace], layout = layout)
        plotly.plot(fig, filename = path)

    def plot_novel_kmer_overlap_percentage(self, tracks):
        path = os.path.join(self.get_current_job_directory(), 'plot_novel_kmer_overlap_percentage.html')
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
        layout = graph_objs.Layout(title = 'Novel kmer Overlap')
        fig = graph_objs.Figure(data = [trace], layout = layout)
        plotly.plot(fig, filename = path)


# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference")
    parser.add_argument("--bed")
    args = parser.parse_args()
    # 
    config.configure(reference_genome = args.reference, bed_file = args.bed)
    #
    novel = HighCoverageNovelJob(job_name = 'novel_', previous_job_name = 'break_point_')
    overlap = NovelKmerOverlapJob(job_name = 'overlap_', previous_job_name = 'novel_')
    #
    novel.execute()