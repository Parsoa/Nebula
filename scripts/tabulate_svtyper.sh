#!/bin/bash
################################################################################
awk 'BEGIN {FS = "\t"} { if ($NF == "00" || $NF == "0/0" || $NF == "./."  || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 00.bed
awk 'BEGIN {FS = "\t"} { if ($NF == "11" || $NF == "1/1" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 11.bed
awk 'BEGIN {FS = "\t"} { if ($NF == "10" || $NF == "0/1" || $NF == "1/0"  || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 10.bed
cat 11.bed > present.bed
cat 10.bed | grep -v CHROM >> present.bed
