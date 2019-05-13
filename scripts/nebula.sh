#!/bin/bash
cd /share/hormozdiarilab/Codes/NebulousSerendipity
cd src/python
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
python -m nebula "$@"
