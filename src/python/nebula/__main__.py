from __future__ import print_function

from nebula import (
    config,
    map_reduce,
)

# import jobs
from cgc import *
from misc import *
#from depth import *
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

def load_tracks():
    job = TrackPreprocessorJob()
    tracks = job.execute()
    config.Configuration.update({'tracks': tracks})

def preprocess():

    def supply_inner_kmers():
        job = ExtractInnerKmersJob()
        job.execute()
        job = ExtractLociIndicatorKmersJob()
        job.execute()
        job = FilterLociIndicatorKmersJob()
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

    load_tracks()
    supply_inner_kmers()
    supply_junction_kmers()
    job = MixKmersJob()
    job.execute()

# ============================================================================================================================ #

def cgc_genotype():
    c = config.Configuration()
    load_tracks()
    job = CgcCounterJob(resume_from_reduce = c.reduce)
    stats = job.execute()
    config.Configuration.update(stats)
    job = CgcIntegerProgrammingJob()
    job.execute()

def genotype():
    c = config.Configuration()
    if c.cgc:
        cgc_genotype()
    else:
        load_tracks()
        job = CgcCounterJob(resume_from_reduce = c.reduce)
        if c.reduce:
            stats = {'coverage': 40, 'std': 8}
        else:
            stats = job.execute()
        config.Configuration.update(stats)
        #job = CgcCoverageCorrectingIntegerProgrammingJob()
        #job.execute()
        #exit()
        job = CgcIntegerProgrammingJob()
        job.execute()
        exit()
        #job = CgcInnerKmersIntegerProgrammingJob()
        #job.execute()
        job = CgcJunctionKmersIntegerProgrammingJob()
        job.execute()

def cluster():
    c = config.Configuration()
    load_tracks()
    job = CgcClusteringJob(begin = 1000, end = 1005)
    job = UnifiedGenotypingJob(begin = 1000, end = 1050, genotyping_batch = 0)
    #job = UnifiedGenotypingOrchestrator()
    job.execute()

def export():
    c = config.Configuration()
    load_tracks()
    job = cgc.ExportGenotypingKmersJob()
    job.execute()

def simulate():
    c = config.Configuration()
    if c.seed == 0:
        print(red('Argument error. Must provide --seed'))
    load_tracks()
    job = Simulation()
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
        if c.command == 'simulate':
            simulate()
        if c.command == 'cluster':
            cluster()
        if c.command == 'export':
            export()
