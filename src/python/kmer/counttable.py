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
)

import khmer
import colorama

print('importing couttable.py')

class DummyCountTable(object):

    def get_kmer_counts(kmer):
        print('dummy')
        return [40]

# ============================================================================================================================ #
# Counttable Import/Export/Creation
# ============================================================================================================================ #

@commons.measure_time
def export_counttable():
    c = config.Configuration()
    # 
    cache = c.counttable + '.ct'
    print(colorama.Fore.BLUE + 'searching for cached counttable ', cache)
    if os.path.isfile(cache):
        print(colorama.Fore.BLUE + 'found at ', cache)
        return
    #
    print(colorama.Fore.BLUE + 'not found, generating counttable...')
    counttable, nkmers = count_kmers_from_file(c.counttable)
    counttable.save(cache)
    print(colorama.Fore.BLUE + 'done')
    #
    print(colorama.Fore.BLUE + 'counttable cached\n', 'kmers: ', nkmers,\
        '\nsize: ', os.stat(cache).st_size)
    return

@commons.measure_time
def import_counttable():
    c = config.Configuration()
    print(colorama.Fore.MAGENTA + 'importing counttable for ', c.counttable)
    cache = c.counttable + '.ct'
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
