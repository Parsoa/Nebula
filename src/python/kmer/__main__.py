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

# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.configure()
    fork.execute()
