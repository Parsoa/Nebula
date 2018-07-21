#!/bin/bash
K=$1
K=$(echo $K | sed 's/C/Z/g')
K=$(echo $K | sed 's/G/C/g')
K=$(echo $K | sed 's/Z/G/g')
K=$(echo $K | sed 's/A/Z/g')
K=$(echo $K | sed 's/T/A/g')
K=$(echo $K | sed 's/Z/T/g')
echo $K | sed 's/./&\n/g' | sed -ne $'x;H;${x;s/\\n//g;p;}'
