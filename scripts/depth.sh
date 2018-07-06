#!/bin/bash
P=$(pwd)
echo $P
SIM=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--simulation " $8 } else { print "" } }')
echo $SIM
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $7 } else { print "hg38.exons.filtered.bed" } }')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
GEN=$(echo $P | awk -F/ '{ print $8 }')
EXN=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "hg38.chr1.exons.filtered.bed" } else { print "hg38.exons.filtered.bed" } }')
echo $BED
echo $REF
echo $GEN
echo $EXN
JLY=$(echo $SIM | awk -v BED=$BED -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/control.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo $JLY
python -m kmer.genotyping --job DepthOfCoverageEstimationJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --exons /share/hormozdiarilab/Codes/NebulousSerendipity/data/$EXN --threads 48 --jellyfish $JLY --fastq $GEN --reference $REF $SIM "$@"
