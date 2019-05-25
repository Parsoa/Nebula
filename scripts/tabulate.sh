#!/bin/bash
################################################################################
awk 'BEGIN {FS = "\t"} { if ($4 ~ /00/)    { print $0 } }' $1 > 00.bed
awk 'BEGIN {FS = "\t"} { if ($4 ~ /11|10/) { print $0 } }' $1 > present.bed
awk 'BEGIN {FS = "\t"} { if ($4 ~ /11/)    { print $0 } }' $1 > 11.bed
awk 'BEGIN {FS = "\t"} { if ($4 ~ /10/)    { print $0 } }' $1 > 10.bed
exit
bedtools intersect -header -a $DATA$BED -b $DATA$TARGET_BED -wa -f 1.0 -F 1.0 -v | sort -u >> absent.bed
#
bedtools intersect -header -a 11.bed -b "homozygous.bed" -wa -f 1.0 -F 1.0 | sort -u > 11_as_11.bed
bedtools intersect -header -a 10.bed -b "homozygous.bed" -wa -f 1.0 -F 1.0 | sort -u > 11_as_10.bed
bedtools intersect -header -a 00.bed -b "homozygous.bed" -wa -f 1.0 -F 1.0 | sort -u > 11_as_00.bed
bedtools intersect -header -a 11.bed -b "heterozygous.bed" -wa -f 1.0 -F 1.0 | sort -u > 10_as_11.bed
bedtools intersect -header -a 10.bed -b "heterozygous.bed" -wa -f 1.0 -F 1.0 | sort -u > 10_as_10.bed
bedtools intersect -header -a 00.bed -b "heterozygous.bed" -wa -f 1.0 -F 1.0 | sort -u > 10_as_00.bed
bedtools intersect -header -a 11.bed -b "absent.bed" -wa -f 1.0 -F 1.0 | sort -u > 00_as_11.bed
bedtools intersect -header -a 10.bed -b "absent.bed" -wa -f 1.0 -F 1.0 | sort -u > 00_as_10.bed
bedtools intersect -header -a 00.bed -b "absent.bed" -wa -f 1.0 -F 1.0 | sort -u > 00_as_00.bed
#
bedtools intersect -header -a 00.bed -b $DATA$TARGET_BED -f 1.0 -F 1.0 > false_negative.bed
sort -n -k2 -o false_negative.bed false_negative.bed

bedtools intersect -header -a 00.bed -b $DATA$TARGET_BED -f 1.0 -F 1.0 -v > true_negative.bed
sort -n -k2 -o true_negative.bed true_negative.bed

bedtools intersect -header -a present.bed -b $DATA$TARGET_BED -f 1.0 -F 1.0 > true_positive.bed
sort -n -k2 -o true_positive.bed true_positive.bed

bedtools intersect -header -a present.bed -b $DATA$TARGET_BED -f 1.0 -F 1.0 -v > false_positive.bed
sort -n -k2 -o false_positive.bed false_positive.bed
wc -l *.bed
