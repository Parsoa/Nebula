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

def preprocess():

    def supply_inner_kmers():
        job = reduction.ExtractInnerKmersJob()
        tracks = job.execute()
        job = reduction.ExtractLociIndicatorKmersJob()
        tracks = job.execute()
        job = reduction.FilterLociIndicatorKmersJob(tracks = tracks)
        job.execute()

    def supply_junction_kmers():
        tracks = {}
        job = junction.ExtractJunctionKmersJob(resume_from_reduce = False)
        job.execute()
        job = junction.JunctionKmersScoringJob(resume_from_reduce = False)
        tracks = job.execute()
        job = junction.FilterJunctionKmersJob(tracks = tracks)
        job.execute()

    load_tracks()
    #supply_inner_kmers()
    #supply_junction_kmers()
    job = preprocessor.MixKmersJob()
    job.execute()

# ============================================================================================================================ #

def cgc_genotype():
    c = config.Configuration()
    load_tracks()
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
        if c.reduce:
            tracks, stats = job.execute()
        else:
            tracks, stats = job.execute()
        config.Configuration.update(stats)
        job = cgc.CgcIntegerProgrammingJob(tracks = tracks)
        job.execute()
        job = cgc.CgcInnerKmersIntegerProgrammingJob(tracks = tracks)
        job.execute()
        job = cgc.CgcJunctionKmersIntegerProgrammingJob(tracks = tracks)
        job.execute()
        if c.select:
            job = cgc.ExportGenotypingKmersJob()
            job.execute()

def verify():
    c = config.Configuration()
    load_tracks()
    job = pacbio.PacBioVerificationJob() 
    job.execute()

def export():
    c = config.Configuration()
    load_tracks()
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
    if c.command == 'export':
        export()
    if c.command == 'verify':
        verify()
    if c.command == 'cluster':
        cluster()
    if c.command == 'simulate':
        simulate()
    if c.command == 'genotype':
        genotype()
    if c.command == 'preprocess':
        preprocess()
