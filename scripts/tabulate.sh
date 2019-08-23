#!/bin/bash
################################################################################
awk 'BEGIN {FS = "\t"} { if ($4 == "00" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 00.bed
awk 'BEGIN {FS = "\t"} { if ($4 != "00" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > present.bed
awk 'BEGIN {FS = "\t"} { if ($4 == "11" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 11.bed
awk 'BEGIN {FS = "\t"} { if ($4 == "10" || $1 == "CHROM" || $1 == "#CHROM") { print $0 } }' $1 > 10.bed
