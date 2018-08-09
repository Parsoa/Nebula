#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source P.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
JLY=$(echo $SIM | awk -v BED=$BED -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo $JLY
echo $@
python -m kmer.programming --job CountInnerKmersJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --fastq "$P/../Simulation/test.fq" --ksize $KSZ $SIM "$@"
