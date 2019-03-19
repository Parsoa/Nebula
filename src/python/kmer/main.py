from __future__ import print_function

from kmer import (
    config,
    map_reduce,
)

from kmer.alu import *
from kmer.misc import *
from kmer.depth import *
from kmer.gapped import *
from kmer.dynamic import *
from kmer.reduction import *
from kmer.simulator import *
from kmer.clustering import *
from kmer.production import *
from kmer.programming import *

# ============================================================================================================================ #

def preprocess():
    job = programming.ExtractInnerKmersJob() 
    job.execute()
    job = reduction.ExtractLociIndicatorKmersJob()
    job.execute()
    job = reduction.FilterLociIndicatorKmersJob()
    job.execute()

# ============================================================================================================================ #

def genotype():
    job = programming.ExtractInnerKmersJob() 
    job.execute()
    job = reduction.ExtractLociIndicatorKmersJob()
    job.execute()
    job = reduction.FilterLociIndicatorKmersJob()
    job.execute()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.reduce)
    # continue
    # if c.command == 'simulate':
    #    Simulation().launch(resume_from_reduce = c.reduce)
