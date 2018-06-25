#!/bin/bash
P=$(pwd)
echo $P
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
BED=$(echo $P | awk -F/ '{print $7}')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
echo $BED
echo $REF
RND=$(echo $BED | awk '{ if ($0 ~ /.*Random.*/) { print "--random" } else { print "" } }')
echo $RND
python -m kmer.simulator --reference $REF --coverage 30 --seed 165784623 --heterozygous $RND --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED
deactivate
