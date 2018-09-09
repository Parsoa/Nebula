#!/bin/bash
source P.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
RND=$(echo $BED | awk '{ if ($0 ~ /.*Random.*/) { print "--random" } else { print "" } }')
python -m kmer.simulator --job $JOB --reference $REF --seed 165784623 --heterozygous $RND --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --jellyfish $JLY $RJF --threads 16 $SIM --chrom $CHR $DESCP $DEPTH "$@"
deactivate
