#!/bin/bash
if [ -z "$1" ]; then
    wc -l 00*.bed 10*.bed 11*.bed
    exit
fi
################################################################################
grep -E "CHROM|0/0" $1 > 00.bed
#awk 'BEGIN {FS = "\t"} { if ($4 == "0/0" || $NF == "0/0" || $4 == "./."  || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 00.bed
grep -E "CHROM|1/1" $1 > 11.bed
#awk 'BEGIN {FS = "\t"} { if ($4 == "1/1" || $NF == "1/1" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 11.bed
grep -E "CHROM|1/0|0/1" $1 > 10.bed
#awk 'BEGIN {FS = "\t"} { if ($4 == "1/0" || $NF == "0/1" || $NF == "1/0" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 10.bed
cat 11.bed > present.bed
cat 10.bed | grep -v CHROM >> present.bed
