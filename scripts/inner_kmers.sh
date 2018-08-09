#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source P.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
JLY=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/control.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo $JLY
source RefJellyfish.sh
BED=$(echo $SIM | awk -v P=$P -v BED=$BED '{ if ($0 ~ /.*kirsimulation.*/) { print P "/../Simulation/all.bed" } else { print "/share/hormozdiarilab/Codes/NebulousSerendipity/data/" BED } }')
echo $BED
python -m kmer.programming --job ExtractInnerKmersJob --bed $BED --threads 48 --reference $REF --jellyfish $JLY $RJF --readlength 100 --insertsize 500 --ksize $KSZ $SIM "$@"
