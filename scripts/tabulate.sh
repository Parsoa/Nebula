#!/bin/bash
################################################################################
awk 'BEGIN {FS = "\t"} { if ($4 == "00")    { print $0 } }' $1 > 00.bed
awk 'BEGIN {FS = "\t"} { if ($4 == "11" || $4 == "10") { print $0 } }' $1 > present.bed
awk 'BEGIN {FS = "\t"} { if ($4 == "11")    { print $0 } }' $1 > 11.bed
awk 'BEGIN {FS = "\t"} { if ($4 == "10")    { print $0 } }' $1 > 10.bed
