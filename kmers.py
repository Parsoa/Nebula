import argparse

def init():
    parse_args()
    exit()

def parse_args():
    parser = argparse.ArgumentParser(add_help = False)
    parser.add_argument("--kmers")
    parser.add_argument("--kmerss")
    parser.parse_args()

if __name__ == '__main__':
    init()
