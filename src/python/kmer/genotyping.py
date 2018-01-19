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
    statistics,
    count_server,
)

# /share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio/HG00732.fq
# /share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio/HG00731.fq
# /share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio/HG00513.fq
# /share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio/HG00512.fq

print('importing genotyping.py')
# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job to produce a kmer signature for each break point of a deletion
# We need to have generated a counttable for the genome we are interested in beforehand
# How to get this far?
# Step 1: get the break points all of which's kmers are found in CHM1, these are the possible candidates for the structural
# variation's boundaries -> BreakPointJob
# Step 2: Find a set of novel kmers for each break point that can be used to indentify it. khmer never underestimates counts so
# if a kmer comes with a count of zero in reference genome, we can be sure that it is really novel -> NovelKmerJob
# Step 3: khmer may report oversetimated counts for these break points so we need to count them exactly again. This is necessary
# for a reliable likelihood model -> CountKmersExactJob
# Step 4: With exact kmer counts available, we can find the most likely break points for each event in our library -> MostLikelyBreakPointsJob
# Step 5: Given a sample genome, try to genotype the structural variations using the likelihood model and signatures gathered
# above -> GenotypingJob (this one)
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenotypingJob(map_reduce.Job):

    # ============================================================================================================================ #
# Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch():
        job = GenotypingJob(job_name = 'Genotyping_', previous_job_name = 'MostLikelyBreakPoints_')
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    # needs to know the most likely breakpoint and its kmers, MostLikelyBreakPointsJob includes that information
    def load_inputs(self):
        self.events = {}
        self.tracks = {}
        path = os.path.join(self.previous_job_directory(), 'merge.json')
        with open(path, 'r') as json_file:
            self.tracks = json.load(path)
            for track in self.tracks:
                self.events[track] = {}
                for break_point in self.tracks[track]['most_likely']:
                    self.events[track][break_point] = self.tracks[track]['break_points'][break_point]
       for index in range(0, self.num_threads):
           self.batch[index] = {}
       index = 0
       for event in self.events:
           self.bathc[index][event] = event
           index = index + 1
           if index == self.num_threads:
               index = 0

    def transform(self, track, track_name):
        c = config.Configuration()
        likelihood = {}
        distribution = {
            (1, 1): statistics.NormalDistribution(mean = c.coverage, c.std = 5),
            (1, 0): statistics.NormalDistribution(mean = c.coverage / 2, c.std = 5),
            (0, 0): statistics.ErrorDistribution(p = 1.0 / 1000),
        }
        # all kmers remaining in the track are to be considered part of the signature, no matter the number of breakpoints
        for break_point in track:
            likelihood[break_point]
            for zyg in distribution:
                likelihood[break_point][zyg] = 0
            for kmer in track[break_point]['novel_kmers']:
                for zyg in distribution:
                    likelihood[break_point][zyg] += distribution[zyg].log_pmf(track[break_point]['novel_kmers'][kmer]['actual_count'])
i       # now find the maximum, for each zygosity find the maximum value, then compare
        m = {}
        for zyg in distribution: 
            for break_point in break_points:
                 if not zyg in m:
                     m[zyg] = lieklihood[break_point][zyg]
                 m[zyg] = max((m[zyg][0], m[zyg][1]), (likelihood[break_point][zyg], break_point), key = operator.itemgetter(0))
        choice = max(m.items(), key = operator.itemgetter(1))[0]
        return {
            'prediction': choice
        }

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    #
    GenotypingJob.launch()
