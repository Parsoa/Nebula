#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source P.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source CovStd.sh
JLY=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo $JLY
source RefJellyfish.sh
echo "$@"
python -m kmer.programming --job IntegerProgrammingJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --jellyfish $JLY $RJF --coverage $COV --fastq $GEN --ksize $KSZ $SIM "$@"
