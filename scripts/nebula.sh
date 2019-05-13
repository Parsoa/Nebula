#!/bin/bash
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
python -m nebula "$@"
