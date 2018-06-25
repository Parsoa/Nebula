import os
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
            path = os.path.join(sys.argv[1], 'most_likely_' + tokens[0] + '_' + tokens[1] + '_' + tokens[2] + '.json')
            n = 0
            with open(path, 'r') as json_file:
                break_points = json.load(json_file)
                for break_point in break_points:
                    start = break_point[1 : break_point.find(',')]
                    end = break_point[break_point.find(',') + 1 : -1]
                    print(start, end)
                    if start != end or break_points[brea]['boundary'] != break_points['(1,1)']['boundary']:
                        bed_file.write(tokens[0] + '\t' + tokens[1] + '\t' + tokens[2] + '\t' + break_point + '\n')
