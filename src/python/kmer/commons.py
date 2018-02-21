import os
import time
import functools
import linecache
import tracemalloc

import colorama

# ============================================================================================================================ #
# Constants
# ============================================================================================================================ #

class bcolors:
    COUNTTABLE = '\033[95m'
    BOUNDARIES = '\033[92m'

# ============================================================================================================================ #
# Decorators
# ============================================================================================================================ #

def measure_time(f):
    def wrapper(*argv):
        start = time.clock()
        result = f(*argv)
        end = time.clock()
        print(colorama.Fore.YELLOW + 'took ', '{:10.9f}'.format(end - start))
        return result
    return wrapper

def display_top(snapshot, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

# ============================================================================================================================ #
# STDIO Wrappers/Helpers
# ============================================================================================================================ #

def colorize(sep = ' ', *vargs):
    return ''.join(functools.reduce(lambda x, y: x + str(y) + sep, vargs))

def red(*args):
    return colorize(colorama.Fore.RED, *args)

def cyan(*args):
    return colorize(colorama.Fore.CYAN, *args)

def blue(*args):
    return colorize(colorama.Fore.BLUE, *args)

def white(*args):
    return colorize(colorama.Fore.WHITE, *args)

def green(*args):
    return colorize(colorama.Fore.GREEN, *args)

def magenta(*args):
    return colorize(colorama.Fore.MAGENTA, *args)

def pretty_print(*args):
    def inner(*vargs):
        return ''.join(functools.reduce(lambda x, y: x + str(y) + ' ', vargs))
    print(inner(colorama.Fore.WHITE, *args))

def json_print(d):
    print(json.dumps(d, sort_keys = True, indent = 4, separators = (',', ': ')))

