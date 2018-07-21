#!/bin/bash
P=$(pwd)
echo $P
SIM=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--simulation " $8 } else { print "" } }')
echo $SIM
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
BED=$(echo $P | awk -F/ '{print $7}')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
CHR=$(echo $BED | awk -F. '{print $3}')
echo $BED
echo $REF
echo $CHR
RND=$(echo $BED | awk '{ if ($0 ~ /.*Random.*/) { print "--random" } else { print "" } }')
echo $RND
python -m kmer.simulator --job Simulation --reference $REF --seed 165784623 --heterozygous $RND --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --jellyfish /share/hormozdiarilab/Experiments/Jellyfish/hg38/mer_counts.jf --threads 8 $SIM --chrom $CHR "$@"
deactivate
