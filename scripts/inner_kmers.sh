#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
P=$(pwd)
echo $P
SIM=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--simulation " $8 } else { print "" } }')
echo $SIM
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
BED=$(echo $P | awk -F/ '{print $7}')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
GEN=$(echo $BED | awk -F. '{ print $1 }')
echo $BED
echo $REF
echo $GEN
DEP=$(echo $P | awk -v P=$P '{ if ($0 ~ /.*simulation.*/) { print P "/../DepthOfCoverageEstimationJob/stats.json"} else { print P "/../../DepthOfCoverageEstimationJob/stats.json" } }')
JLY=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/control.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo $JLY
python -m kmer.programming --job ExtractInnerKmersJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --reference $REF --jellyfish $JLY /share/hormozdiarilab/Experiments/Jellyfish/$REF/mer_counts.jf --readlength 100 --insertsize 500 $SIM "$@"
