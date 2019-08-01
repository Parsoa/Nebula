#!/bin/bash
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
echo Launching Nebula...
python -m nebula "$@"
