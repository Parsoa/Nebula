#!/bin/bash
P=$(pwd)
echo $P
SIM=$(echo $P | awk '{ if ($0 ~ /.*mulation.*/) { print "--simulation" } else { print "" } }')
echo $SIM
BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $7 } else { print $10 } }')
echo $BED
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/) { print "hg38" } else { print "hg19" } }')
echo $REF
GEN=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "SIM" } else { print $8 } }')
echo $GEN
BED_GENOME=$(echo $BED | awk -F "." '{print $1}')
echo $BED_GENOME
EVENT=$(pwd | awk -F "/" '{ if ($0 ~ /'.*DEL.*'/) { print "DEL" } else { print "INV" } }' | xargs)
echo $EVENT
DATA="/share/hormozdiarilab/Codes/NebulousSerendipity/data/"
TARGET=$GEN.$REF.$EVENT.bed
echo $TARGET
################################################################################
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/11.bed -b ./11.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_11_as_11.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/11.bed -b ./10.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_11_as_10.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/11.bed -b ./00.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_11_as_00.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/10.bed -b ./11.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_10_as_11.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/10.bed -b ./10.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_10_as_10.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/10.bed -b ./00.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_10_as_00.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/00.bed -b ./11.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_00_as_11.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/00.bed -b ./10.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_00_as_10.bed
bedtools intersect -a ../GappedKmersIntegerProgrammingJob/00.bed -b ./00.bed -wa -f 1.0 -F 1.0 | sort -u > gapped_00_as_00.bed
################################################################################
