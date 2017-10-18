import os
import pwd

from sys import platform

from kmer import (
    bed,
    fork,
    sets,
    config,
    reference,
)

import khmer
import colorama
import pybedtools

# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.configure()
    fork.execute()
