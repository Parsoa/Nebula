export BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*genotyping.*/) { print $10 } else { print $7 } }')
echo BED $BED
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/) { print "hg38" } else { print "hg19" } }')
echo REF $REF
GEN=$(echo $BED | awk -F. '{ print $1 }')
echo GEN $GEN
CHR=$(echo $BED | awk -F. '{ print $3 }')
echo CHR $CHR
