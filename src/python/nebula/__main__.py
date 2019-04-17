from __future__ import print_function

from nebula import (
    config,
    map_reduce,
)

# import jobs
from alu import *
from cgc import *
from misc import *
from depth import *
from gapped import *
from junction import *
from reduction import *
from simulator import *
from clustering import *
from production import *
from programming import *
from preprocessor import *

# import helpers
from debug import *
from nebulas import *
from commons import *
from chromosomes import *

print = pretty_print

# ============================================================================================================================ #

def preprocess():

    def supply_inner_kmers():
        job = ExtractInnerKmersJob()
        job.execute()
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
        #job = junction.ExtractJunctionKmersJob(resume_from_reduce = False)
        #job.execute()
        job = UniqueJunctionKmersJob()
        job.execute()
        job = JunctionKmersScoringJob()
        job.execute()
        job = FilterJunctionKmersJob()
        job.execute()

    job = TrackPreprocessorJob()
    tracks = job.execute()
    config.Configuration.update({'tracks': tracks})
    #extract_whole_genome()
    #supply_inner_kmers()
    #supply_gapped_kmers()
    supply_junction_kmers()
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
