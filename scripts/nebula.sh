#!/bin/bash
WORKDIR=$PWD
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv3/bin/activate
cd src/python
#export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
if [ "$(hostname)" == "rum" ]; then
    export PYTHONPATH=$PYTHONPATH:/share/hormozdiarilab/Codes/NebulousSerendipity/Jellyfish/install/lib/python3.5/site-packages
    python -m nebula "$@"
else
    python -m nebula "$@" --workdir $WORKDIR
fi
