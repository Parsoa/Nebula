import time

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

# ============================================================================================================================ #
# STDIO Wrappers/Helpers
# ============================================================================================================================ #

def identitiy(*vargs, sep = ' '):
    return ''.join(functools.reduce(lambda x, y: x + str(y) + sep, vargs))

def white(*args):
    return identitiy(colorama.Fore.WHITE, *args)

def green(*args):
    return identitiy(colorama.Fore.GREEN, *args)

def red(*args):
    return identitiy(colorama.Fore.RED, *args)

def cyan(*args):
    return identitiy(colorama.Fore.CYAN, *args)

def blue(*args):
    return identitiy(colorama.Fore.BLUE, *args)

def magenta(*args):
    return identitiy(colorama.Fore.MAGENTA, *args)

def pretty_print(*args):
    def inner(*vargs):
        return ''.join(functools.reduce(lambda x, y: x + str(y) + ' ', vargs))
    print(inner(colorama.Fore.WHITE, *args))

def json_print(d):
    print(json.dumps(d, sort_keys = True, indent = 4, separators = (',', ': ')))

