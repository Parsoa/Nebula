from __future__ import print_function

from nebula import (
    cgc,
    config,
    pacbio,
    junction,
    reduction,
    simulator,
    map_reduce,
    programming,
    preprocessor,
)

from nebula.debug import *
from nebula.kmers import *
from nebula.logger import *
from nebula.chromosomes import *

# ============================================================================================================================ #

def load_tracks(filter_overlap = True):
    job = preprocessor.TrackPreprocessorJob(filter_overlap = filter_overlap)
    tracks = job.execute()
    config.Configuration.update({'tracks': tracks})

def misc():
    job = cgc.PcaClusteringJob()
    job.execute()

def unify():
    job = preprocessor.EventSetUnificationJob()
    job.execute()

def preprocess():

    def supply_inner_kmers():
        c = config.Configuration()
        tracks = {}
        job = reduction.ExtractInnerKmersJob()
        tracks = job.execute()
        job = reduction.ScoreInnerKmersJob()
        tracks = job.execute()
        job = reduction.FilterInnerKmersJob(tracks = tracks)
        job.execute()

    def supply_junction_kmers():
        c = config.Configuration()
        tracks = {}
        job = junction.ExtractJunctionKmersJob()
        job.execute()
        if c.cpp:
            return
        job = junction.ScoreJunctionKmersJob()
        tracks = job.execute()
        job = junction.FilterJunctionKmersJob(tracks = tracks)
        job.execute()

    load_tracks()
    #supply_inner_kmers()
    supply_junction_kmers()
    preprocessor.MixKmersJob().execute()

# ============================================================================================================================ #

def cgc_genotype():
    c = config.Configuration()
    load_tracks(filter_overlap = False)
    job = cgc.CgcCounterJob(resume_from_reduce = c.reduce)
    tracks, stats = job.execute()
    config.Configuration.update(stats)
    job = cgc.CgcIntegerProgrammingJob(tracks = tracks)
    job.execute()

def genotype():
    c = config.Configuration()
    if c.cgc:
        cgc_genotype()
    else:
        load_tracks(filter_overlap = False)
        job = cgc.CgcCounterJob(resume_from_reduce = c.reduce)
        #job.execute()
        job = cgc.CgcIntegerProgrammingJob()
        job.execute()
        if c.select:
            job = cgc.ExportGenotypingKmersJob()
            job.execute()

def simulate():
    c = config.Configuration()
    if c.seed == 0:
        user_print_warning('Default randomness seed of 0 is used.')
    load_tracks()
    job = simulator.Simulation()
    job.execute()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

def print_banner():
    print('Nebula, ultra-efficient mapping-free structural variation genotyper.')

if __name__ == '__main__':
    config.init()
    print_banner()
    c = config.Configuration()
    if c.command == 'gc':
        gc()
    if c.command == 'misc':
        misc()
    if c.command == 'unify':
        unify()
    if c.command == 'simulate':
        simulate()
    if c.command == 'genotype':
        genotype()
    if c.command == 'preprocess':
        preprocess()
