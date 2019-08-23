from __future__ import print_function

from nebula import (
    config,
    map_reduce,
)

from cgc import *
from misc import *
from junction import *
from reduction import *
from simulator import *
from clustering import *
from programming import *
from preprocessor import *

from debug import *
from kmers import *
from logger import *
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
        tracks = job.execute()
        job = ExtractLociIndicatorKmersJob()
        tracks = job.execute()
        job = FilterLociIndicatorKmersJob(tracks = tracks)
        job.execute()

    def supply_junction_kmers():
        tracks = {}
        #job = junction.ExtractJunctionKmersJob(resume_from_reduce = False)
        #job.execute()
        job = JunctionKmersScoringJob(resume_from_reduce = True)
        tracks = job.execute()
        job = FilterJunctionKmersJob(tracks = tracks)
        job.execute()

    load_tracks()
    #supply_inner_kmers()
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
            stats = {'coverage': 23, 'std': 7}
            stats = {'coverage': 50, 'std': 17}
            tracks, stats = job.execute()
        else:
            tracks, stats = job.execute()
        config.Configuration.update(stats)
        job = CgcIntegerProgrammingJob(tracks = tracks)
        job.execute()
        exit()
        job = CgcInnerKmersIntegerProgrammingJob(tracks = tracks)
        job.execute()
        job = CgcJunctionKmersIntegerProgrammingJob(tracks = tracks)
        job.execute()
        job = ExportGenotypingKmersJob()
        job.execute()

def cluster():
    c = config.Configuration()
    load_tracks()
    #job = CgcClusteringJob(begin = 1000, end = 1005)
    #job = UnifiedGenotypingJob(begin = 1000, end = 1050, genotyping_batch = 0)
    #job = UnifiedGenotypingOrchestrator()
    job = UnifiedGenotypingJob(begin = 0, end = 5, genotyping_batch = 'for-august-presentation')
    job.execute()
    exit()
    for i in range(1, 10):
        job = UnifiedGenotypingJob(begin = i * 5, end = (i + 1) * 5, genotyping_batch = i)
        job.execute()
        exit()

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

def print_banner():
    print('Nebula, ultra-efficient mapping-free structural variation genotyper')

if __name__ == '__main__':
    config.init()
    print_banner()
    c = config.Configuration()
    if c.command == 'export':
        export()
    if c.command == 'cluster':
        cluster()
    if c.command == 'simulate':
        simulate()
    if c.command == 'genotype':
        genotype()
    if c.command == 'preprocess':
        preprocess()
