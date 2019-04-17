from __future__ import print_function

from kmer import (
    config,
    map_reduce,
)

# import jobs
from alu import *
from cgc import *
from misc import *
from depth import *
from gapped import *
from dynamic import *
from reduction import *
from simulator import *
from clustering import *
from production import *
from programming import *

# import helpers
from debug import *
from kmers import *
from commons import *
from chromosomes import *

print = pretty_print

# ============================================================================================================================ #

def preprocess():

    def supply_inner_kmers():
        job = programming.ExtractInnerKmersJob()
        job.execute()
        job = reduction.ExtractLociIndicatorKmersJob()
        job.execute()
        job = reduction.FilterLociIndicatorKmersJob()
        job.execute()

    def supply_gapped_kmers():
        job = gapped.ExtractGappedKmersJob()
        job.execute()
        job = gapped.UniqueGappedKmersJob()
        job.execute()
        job = gapped.UniqueGappedKmersScoringJob()
        job.execute()

    def supply_junction_kmers():
        job = dynamic.ExtractJunctionKmersJob()
        job.execute()
        job = dynamic.UniqueJunctionKmersJob()
        job.execute()
        job = dynamic.JunctionKmersScoringJob()
        job.execute()
        job = dynamic.FilterJunctionKmersJob()
        job.execute()

    #extract_whole_genome()
    #supply_inner_kmers()
    #supply_gapped_kmers()
    #supply_junction_kmers()
    job = production.MixKmersJob()
    job.execute()

# ============================================================================================================================ #

def genotype():
    c = config.Configuration()
    if c.cgc:
        job = CgcCounterJob(resume_from_reduce = c.reduce)
        tracks, stats = job.execute()
        stats['threads'] = 4
        config.Configuration.update(stats)
        job = CgcIntegerProgrammingJob()
        job.tracks = tracks
        job.execute()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    'Nebula, ultra-efficient mapping-free structural variation genotyper'
    config.init()
    c = config.Configuration()
    if c.command == 'preprocess':
        preprocess()
    if c.command == 'genotype':
        genotype()
    #getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.reduce)
    # continue
    # if c.command == 'simulate':
    #    Simulation().launch(resume_from_reduce = c.reduce)
