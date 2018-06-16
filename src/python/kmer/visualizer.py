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

import plotly.offline as plotly
import plotly.graph_objs as graph_objs

import numpy

# ============================================================================================================================ #
# Plotly helpers
# ============================================================================================================================ #

def histogram(x, name, path):
    data = [graph_objs.Histogram(x = x, xbins = dict(start = min(x), size = 5, end = max(x)))]
    plotly.plot(data, filename = os.path.join(path, 'histogram_' + name + '.html'), auto_open = False)
