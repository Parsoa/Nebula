#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
P=$(pwd)
echo $P
SIM=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--simulation " $8 } else { print "" } }')
echo $SIM
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $7 } else { print $10 } }')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
GEN=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "SIM" } else { print $8 } }')
echo $BED
echo $REF
echo $GEN
JLY=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo $JLY
echo "$@"
python -m kmer.depth --job UniqueKmersDepthOfCoverageEstimationJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --jellyfish $JLY /share/hormozdiarilab/Experiments/Jellyfish/$REF/mer_counts.jf --fastq "$P/../Simulation/test.fq" $SIM "$@"
