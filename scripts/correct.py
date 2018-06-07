import sys
import json

import plotly.offline as plotly
import plotly.graph_objs as go

target = sys.argv[2]
print(target)
with open('./' + target + '.bed') as bed_file:
    lines = bed_file.readlines()[1:]
    with open('./incorrect_' + target + '.bed', 'w') as bed_file:
        for line in lines:
            tokens = line.split()
            path = '/share/hormozdiarilab/Codes/NebulousSerendipity/output/' + sys.argv[1] + '/31/MostLikelyBreakPointsJob/most_likely_' + name + '.json'
            n = 0
            with open(path, 'r') as json_file:
                break_points = json.load(json_file)['break_points']
                for break_point in break_points:
                    # check if the chosen breakpoint is the same as the actual one
                    start = break_point[0 : break_point.find(',')]
                    end = break_point[break_point.find(',') + 1 : -1]
                    if '(0,0)' not in break_points or start != end or break_points[most_likely]['boundary'] != break_points['(0,0)']['boundary']:
                        bed_file.write(line)
