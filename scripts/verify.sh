bedtools intersect -a 00.bed -b $1 -f 1.0 -F 1.0 > false_negative.bed
bedtools intersect -a 00.bed -b $1 -f 1.0 -F 1.0 -v > true_negative.bed
bedtools intersect -a present.bed -b $1 -f 1.0 -F 1.0 > true_positive.bed
bedtools intersect -a present.bed -b $1 -f 1.0 -F 1.0 -v > false_positive.bed
