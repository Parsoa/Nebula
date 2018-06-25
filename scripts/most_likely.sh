#!/bin/bash
P=$(pwd)
echo $P
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
BED=$(echo $P | awk -F/ '{print $7}')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
GEN=$(echo $BED | awk -F. '{ print $1 }')
echo $BED
echo $REF
echo $GEN
DEP=$(echo $P | awk -v P=$P '{ if ($0 ~ /.*simulation.*/) { print P "/../DepthOfCoverageEstimationJob/stats.json"} else { print P "/../../DepthOfCoverageEstimationJob/merge.json" } }')
echo $DEP
COV="$(head -n 10 "$DEP" | sed -n 's/\s*"mean":\s*\([0-9]*\).*/\1/p')"
STD="$(head -n 10 "$DEP" | sed -n 's/\s*"std":\s*\([0-9]*\).*/\1/p')"
echo $COV
echo $STD
SIM=$(echo $P | awk '{ if ($0 ~ /.*mulation.*/) { print "--simulation" } else { print "" } }')
echo $SIM
python -m kmer.break_point --job MostLikelyBreakPointsJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --coverage $COV --std $STD $SIM "$@"
