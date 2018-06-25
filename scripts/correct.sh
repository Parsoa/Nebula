#!/bin/bash
P=$(pwd)
echo $P
SIM=$(echo $P | awk '{ if ($0 ~ /.*mulation.*/) { print "--simulation" } else { print "" } }')
echo $SIM
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd scripts
BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $7 } else { print $10 } }')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/) { print "hg38" } else { print "hg19" } }')
GEN=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "SIM" } else { print $8 } }')
CORRECT=$(readlink -e correct.py)
MST=$(echo $P | awk -v P=$P -v BED=$BED '{ if ($0 ~ /.*mulation.*/) { print P "/../" } else { print P "/../../../../ BED "/31/" } }')
echo $MST
cd $P
exit
python $CORRECT $MST "merge"
