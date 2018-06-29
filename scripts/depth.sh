#!/bin/bash
P=$(pwd)
echo $P
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
SIM=$(echo $P | awk '{ if ($0 ~ /.*mulation.*/) { print "--simulation" } else { print "" } }')
echo $SIM
BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $7 } else { print "hg38.exons.filtered.bed" } }')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
GENOME=$(echo $P | awk -F/ '{ print $8 }')
EXN=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "hg38.chr1.exons.filtered.bed" } else { print "hg38.exons.filtered.bed" } }')
echo $BED
echo $REF
echo $GENOME
echo $EXN
JLY=$(echo $SIM | awk -v BED=$BED -v P=$P -v GENOME=$GENOME '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/control.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GENOME "/mer_counts.jf" } }')
echo $JLY
python -m kmer.genotyping --job DepthOfCoverageEstimationJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --exons /share/hormozdiarilab/Codes/NebulousSerendipity/data/$EXN --threads 48 --jellyfish $JLY --fastq $GENOME --reference $REF $SIM "$@"
