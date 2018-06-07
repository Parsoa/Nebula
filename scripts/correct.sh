#!/bin/bash
PWD=$(pwd)
echo $PWD
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd scripts
CORRECT=$(readlink -e correct.py)
cd $PWD
python $CORRECT `pwd | awk -F "/" '{print $NF}'` "true_positive"
python $CORRECT `pwd | awk -F "/" '{print $NF}'` "true_negative"
python $CORRECT `pwd | awk -F "/" '{print $NF}'` "false_positive"
python $CORRECT `pwd | awk -F "/" '{print $NF}'` "false_negative"
bedtools intersect -a true_positive.bed  -b incorrect_true_positive.bed  -v > correct_true_positive.bed
bedtools intersect -a false_positive.bed -b incorrect_false_positive.bed -v > correct_false_positive.bed
bedtools intersect -a true_negative.bed  -b incorrect_true_negative.bed  -v > correct_true_negative.bed
bedtools intersect -a false_negative.bed -b incorrect_false_negative.bed -v > correct_false_negative.bed
