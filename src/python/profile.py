from kmer import (
    prune,
    config
)

config.init()
prune.CountBedKmersExactJob.launch()
