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

from kmer import (
    config,
)
import pandas as pd
import plotly.offline as plotly
import plotly.graph_objs as graph_objs
#import matplotlib.pyplot as plt
import plotly.figure_factory as ff

import numpy

# ============================================================================================================================ #
# Plotly helpers
# ============================================================================================================================ #

#def matplotlib_histogram(x, name, path):
#    fig, axis = plt.suplots()
#    axis.hist(x)
#    fig.savefig(os.path.join(path, 'histogram_' + name))

def histogram(x, name, path, x_label, y_label, step = 5):
    data = [graph_objs.Histogram(x = x, xbins = dict(start = min(x), size = 1, end = max(x)))]
    layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    figure = graph_objs.Figure(data = data, layout = layout)
    plotly.plot(figure, filename = os.path.join(path, 'histogram_' + name + '.html'), auto_open = False)

def scatter(x, y, name, path, x_label, y_label):
    data = [graph_objs.Scatter(x = x, y = y, mode = 'lines')]
    layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    figure = graph_objs.Figure(data = data, layout = layout)
    plotly.plot(figure, filename = os.path.join(path, 'scatter_' + name + '.html'), auto_open = False)

def bar(x, ys, name, path, x_label, y_label):
    data = []
    for y in ys:
        data.append(graph_objs.Bar(x = x, y = y))
    layout = graph_objs.Layout(title = name, xaxis = dict(title = x_label), yaxis = dict(title = y_label))
    figure = graph_objs.Figure(data = data, layout = layout)
    plotly.plot(figure, filename = os.path.join(path, 'bar_' + name + '.html'), auto_open = False)

def violin(x, y, name, path, x_label, y_label):
    data = pd.DataFrame(dict(LP = y, Zygosity = x))
    fig = ff.create_violin(data, data_header = 'LP', group_header = 'Zygosity', width = 400 * len(set(x)))
    plotly.plot(fig, filename = os.path.join(path, 'violin_' + name + '.html'), auto_open = False)

