export BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*genotyping.*/) { print $9 } else { print $7 } }')
echo BED $BED
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/) { print "hg38" } else { print "hg19" } }')
echo REF $REF
if [[ $P =~ "genotyping" ]]; then
    echo YES
    GEN=$(echo $P | awk -F/ '{ print $8 }')
else
    echo NO
    GEN=$(echo $BED | awk -F. '{ print $1 }')
fi
echo GEN $GEN
CHR=$(echo $BED | awk -F. '{ print $3 }')
echo CHR $CHR
