#!/bin/bash
P=$(pwd)
echo $P
SIM=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--simulation " $8 } else { print "" } }')
echo $SIM
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $7 } else { print "hg38.exons.filtered.bed" } }')
CHR=$(echo $BED | awk -F. -v P=$P '{ if (P ~ /.*mulation.*/) { print $2 } else { print "" } }')
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
GEN=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "SIM" } else { print $8 } }')
EXN=$(echo $P | awk -F/ -v CHR=$CHR '{ if ($0 ~ /.*mulation.*/) { print "hg38." CHR  ".exons.filtered.bed" } else { print "hg38.exons.filtered.bed" } }')
echo $BED
echo $CHR
echo $REF
echo $GEN
echo $EXN
JLY=$(echo $SIM | awk -v BED=$BED -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/control.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo $JLY
FSQ=$(echo $SIM | awk -v BED=$BED -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.fq" } else { print "/share/hormozdiarilab/Codes/NebulousSerendipity/data/" GEN ".Fq" } }')
echo $FSQ
echo $@
python -m kmer.depth --job ExonDepthOfCoverageCountingJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --exons /share/hormozdiarilab/Codes/NebulousSerendipity/data/$EXN --threads 48 --jellyfish $JLY --fastq $FSQ --reference $REF $SIM "$@"
