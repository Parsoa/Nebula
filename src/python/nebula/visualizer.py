from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import random
import argparse

from nebula import (
    config,
    statistics,
)

#import pandas as pd
#import plotly.offline as plotly
#import plotly.graph_objs as graph_objs
#import plotly.figure_factory as ff

# ============================================================================================================================ #
# Plotly helpers
# ============================================================================================================================ #

def histogram(x, name, path, x_label, y_label, step = 1):
    #data = [graph_objs.Histogram(x = x, xbins = dict(start = min(x), size = step, end = max(x) + 1))]
    #layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    #figure = graph_objs.Figure(data = data, layout = layout)
    #plotly.plot(figure, filename = os.path.join(path, 'histogram_' + name + '.html'), auto_open = False)
    pass

def scatter(x, y, name, path, x_label, y_label):
    #data = [graph_objs.Scatter(x = x, y = y, mode = 'lines')]
    #layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    #figure = graph_objs.Figure(data = data, layout = layout)
    #plotly.plot(figure, filename = os.path.join(path, 'scatter_' + name + '.html'), auto_open = False)
    pass

def bar(x, ys, name, path, x_label, y_label):
    #data = []
    #for y in ys:
    #    data.append(graph_objs.Bar(x = x, y = y))
    #layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    #figure = graph_objs.Figure(data = data, layout = layout)
    #plotly.plot(figure, filename = os.path.join(path, 'bar_' + name + '.html'), auto_open = False)
    pass

def violin(x, y, name, path, x_label, y_label):
    #xs = []
    #ys = []
    #m = statistics.mean(y)
    #for index, d in enumerate(y):
    #    if d <= 3 * m or d < 200:
    #        xs.append(x[index])
    #        ys.append(d)
    #data = pd.DataFrame(dict(LP = ys, Zygosity = xs))
    #fig = ff.create_violin(data, data_header = 'LP', group_header = 'Zygosity', width = 200 * len(set(x)))
    #plotly.plot(fig, filename = os.path.join(path, 'violin_' + name + '.html'), auto_open = False)
    pass

