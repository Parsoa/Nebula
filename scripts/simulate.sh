#!/bin/bash
source P.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
RND=$(echo $BED | awk '{ if ($0 ~ /.*Random.*/) { print "--random" } else { print "" } }')
export MEM="60G"
export TME="12:00:00"
source SLURM.sh
#python -m kmer.simulator --job $JOB --reference $REF --seed 165784623 --heterozygous $RND --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --jellyfish /share/hormozdiarilab/Experiments/Jellyfish/hg38/mer_counts.jf --threads 16 $SIM --chrom $CHR "$@" deactivate
