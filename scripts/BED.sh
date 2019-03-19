BED=$(echo $P | awk -F/ '{ if ($0 ~ /.*genotyping.*/) { print $(NF - 1) } else if ($0 ~ /.*simulation.*/) { print $(NF - 3) } else { print $(NF - 1)} }')
echo BED $BED
REF=$(echo $BED | awk -F. '{ if ($0 ~ /.*hg38.*/) { print "/share/hormozdiarilab/Data/ReferenceGenomes/Hg38/GRC38.fasta" } else { print "/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref" } }')
echo REF $REF
if [[ $P =~ "genotyping" ]]; then
    echo YES
    GEN=$(echo $P | awk -F/ '{ print $(NF - 2) }')
else
    echo NO
    GEN=$(echo $BED | awk -F. '{ print $1 }')
fi
echo GEN $GEN
CHR=$(echo $BED | awk -F. '{ print $3 }')
echo CHR $CHR
