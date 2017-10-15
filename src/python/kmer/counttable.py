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
def export_sample_counttable():
    c = config.Configuration()
    # 
    print(colorama.Fore.BLUE + 'searching for cached counttable ...')
    cache = c.fastq_file + '.ct'
    if os.path.isfile(cache):
        print(colorama.Fore.BLUE + 'found at ', cache)
        return
    #
    print(colorama.Fore.BLUE + 'not found, generating counttable ' + cache)
    counttable, nkmers = count_kmers_from_file(c.fastq_file)
    counttable.save(cache)
    print(colorama.Fore.BLUE + 'done')
    #
    print(colorama.Fore.BLUE + 'sample counttable cached\n', 'kmers: ', nkmers,\
        '\nsize: ', os.stat(cache).st_size)
    return

@commons.measure_time
def import_sample_counttable():
    print(colorama.Fore.MAGENTA + 'importing counttable ...')
    c = config.Configuration()
    cache = c.fastq_file + '.ct'
    counttable = khmer.Counttable.load(cache)
    print(colorama.Fore.MAGENTA + 'done')
    return counttable

def count_kmers_from_file(file):
    c = config.Configuration()
    #
    counttable = khmer.Counttable(c.ksize, c.khmer_table_size, c.khmer_num_tables)
    nseqs, nkmers = counttable.consume_seqfile(file)
    #
    return counttable, nkmers
