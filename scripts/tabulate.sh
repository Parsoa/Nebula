#!/bin/bash
################################################################################
awk 'BEGIN {FS = "\t"} { if ($4 == "00" || $4 == "0/0" || $4 == "./."  || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 00.bed
awk 'BEGIN {FS = "\t"} { if ($4 == "11" || $4 == "1/1" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 11.bed
awk 'BEGIN {FS = "\t"} { if ($4 == "10" || $4 == "0/1" || $4 == "1/0"  || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 10.bed
cat 11.bed > present.bed
cat 10.bed | grep -v CHROM >> present.bed
