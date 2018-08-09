export BED=$(echo $P | awk -F/ '{print $7}')
echo $BED
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/){print "hg38"} else {print "hg19"}}')
echo $REF
GEN=$(echo $BED | awk -F. '{ print $1 }')
echo $GEN
CHR=$(echo $BED | awk -F. '{print $3}')
echo $CHR
