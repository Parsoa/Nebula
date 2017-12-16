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

# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job to produce a kmer signature for each break point of a deletion
# Step 1: get the break points all of which's kmers are found in CHM1, these are the possible candidates for the structural
# variation's boundaries -> BreakPointJob
# Step 2: Find a set of novel kmers for each break point that can be used to indentify it. khmer never underestimates counts so
# if a kmer comes with a count of zero in reference genome, we can be sure that it is really novel -> NovelKmerJob
# Step 3: khmer may report oversetimated counts for these break points so we need to count them exactly again. This is necessary
# for a reliable likelihood model -> CountKmersExactJob
# Step 4: With exact kmer counts available, we can find the most likely break points for each event in our library -> MostLikelyBreakPointsJob
# Step 5: Given a sample genome, try to genotype the structural variations using the likelihood model and signatures gathered
# above.
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

    def transform(self, track, track_name):
        c = config.Configuration()
        likelihood = {}
        # TODO: proper value for std?
        distribution = {
            (1, 1): statistics.NormalDistribution(mean = c.coverage, std = 5),
            (1, 0): statistics.NormalDistribution(mean = c.coverage / 2, std = 5),
            (0, 0): statistics.NormalDistribution(mean = c.coverage / 2, std = 5),
        }
        zygosity = [(1, 1), (1, 0)]
        break_points = []
        for kmer in track:
            for break_point in track[kmer]['break_points']:
                if not break_point in likelihood:
                    likelihood[break_point] = {
                        (1, 1): 1,
                        (1, 0): 1,
                    }
                    break_points.append(break_point)
                for zyg in zygosity:
                    likelihood[break_point][zyg] *= distribution[zyg](track[kmer]['actual_count'])
        # TODO: each term should be multiplied by P(zyg | pb) , how to calculate
        output = map(lambda x: likelihood[x][(1, 1)] + likelihood[x](1, 0), break_points)
        return output

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #