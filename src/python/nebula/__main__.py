from __future__ import print_function

from nebula import (
    config,
    map_reduce,
)

# import jobs
from cgc import *
from misc import *
from depth import *
from gapped import *
from junction import *
from reduction import *
from simulator import *
from clustering import *
from programming import *
from preprocessor import *

# import helpers
from debug import *
from kmers import *
from commons import *
from chromosomes import *

print = pretty_print

# ============================================================================================================================ #

def preprocess():

    def supply_inner_kmers():
        job = ExtractInnerKmersJob()
        job.execute()
        return
        job = ExtractLociIndicatorKmersJob()
        job.execute()
        job = FilterLociIndicatorKmersJob()
        job.execute()

    def supply_gapped_kmers():
        job = ExtractGappedKmersJob()
        job.execute()
        job = UniqueGappedKmersJob()
        job.execute()
        job = UniqueGappedKmersScoringJob()
        job.execute()

    def supply_junction_kmers():
        job = junction.ExtractJunctionKmersJob(resume_from_reduce = False)
        job.execute()
        job = UniqueJunctionKmersJob(resume_from_reduce = False)
        job.execute()
        job = JunctionKmersScoringJob()
        job.execute()
        job = FilterJunctionKmersJob()
        job.execute()

    job = TrackPreprocessorJob()
    tracks = job.execute()
    config.Configuration.update({'tracks': tracks})
    #supply_inner_kmers()
    #supply_gapped_kmers()
    supply_junction_kmers()
    job = MixKmersJob()
    job.execute()

# ============================================================================================================================ #

def genotype():
    c = config.Configuration()
    job = TrackPreprocessorJob()
    stats = {}
    tracks = job.execute()
    config.Configuration.update({'tracks': tracks})
    job = CgcCounterJob(resume_from_reduce = c.reduce)
    stats = job.execute()
    config.Configuration.update(stats)
    if c.cgc:
        stats['threads'] = c.threads
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
    print('Nebula, ultra-efficient mapping-free structural variation genotyper')
    config.init()
    c = config.Configuration()
    if c.job:
        getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.reduce)
    else:
        if c.command == 'preprocess':
            preprocess()
        if c.command == 'genotype':
            genotype()
