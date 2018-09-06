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
TARGET=$(echo $P | awk '{if ($0 ~ /.*MixGappedUniqueInnerKmersIntegerProgrammingJob.*/) { print "MixUniqueInnerKmersIntegerProgrammingJob" } if ($0 ~ /.*MixGappedUniqueInnerNoneUniqueInnerKmersIntegerProgrammingJob.*/) { print "MixGappedUniqueInnerKmersIntegerProgrammingJob" } }')
################################################################################
pwd
zyg=("11" "10" "00")
for i in "${zyg[@]}"; do
    for j in "${zyg[@]}"; do
        for k in "${zyg[@]}"; do
            bedtools intersect -a ../$TARGET/"$i"_as_"$j".bed -b ./"$i"_as_"$k".bed -wa -f 1.0 -F 1.0 | sort -u > from_"$i"_as_"$j"_to_"$k".bed
            bedtools intersect -a ../$TARGET/no_signal.bed -b ./"$i"_as_"$k".bed -wa -f 1.0 -F 1.0 | sort -u > none_to_"$i"_as_"$k".bed
        done
    done
done
################################################################################
