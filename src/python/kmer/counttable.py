import io
import os
import pwd
import sys
import json
import time

from . import (
    sets,
    config,
    commons,
    reference,
)

import khmer
import colorama

# ============================================================================================================================ #
# Counttable Import/Export/Creation
# ============================================================================================================================ #

@commons.measure_time
def export_counttable(seq_file):
    c = config.Configuration()
    # 
    cache = seq_file + '.ct'
    print(colorama.Fore.BLUE + 'searching for cached counttable ', cache)
    if os.path.isfile(cache):
        print(colorama.Fore.BLUE + 'found at ', cache)
        return
    #
    print(colorama.Fore.BLUE + 'not found, generating counttable...')
    counttable, nkmers = count_kmers_from_file(seq_file)
    counttable.save(cache)
    print(colorama.Fore.BLUE + 'done')
    #
    print(colorama.Fore.BLUE + 'counttable cached\n', 'kmers: ', nkmers,\
        '\nsize: ', os.stat(cache).st_size)
    return

@commons.measure_time
def import_sample_counttable(seq_file):
    print(colorama.Fore.MAGENTA + 'importing counttable for ', seq_file)
    c = config.Configuration()
    cache = seq_file + '.ct'
    counttable = khmer.Counttable.load(cache)
    print(colorama.Fore.MAGENTA + 'done')
    return counttable

def count_kmers_from_file(seq_file):
    c = config.Configuration()
    #
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counttable.consume_seqfile(seq_file)
    #
    return counttable, nkmers
