#!/usr/bin/python

import os
import sys

header = None

def load_tracks(path):
    tracks = {}
    with open(path) as bed_file:
        for line in bed_file.readlines():
            tokens = line.split()
            tracks[tokens[0] + '_' + tokens[1] + '_' + tokens[2]] = '\t'.join(tokens[:3]) + '\n'
    return tracks

a_tracks = load_tracks(os.path.join(sys.argv[1]))
b_tracks = load_tracks(os.path.join(sys.argv[2]))

tracks = {}

for track in a_tracks:
    tracks[track] = a_tracks[track]

for track in b_tracks:
    if not track in tracks:
        tracks[track] = b_tracks[track]

for track in tracks:
    print tracks[track]
