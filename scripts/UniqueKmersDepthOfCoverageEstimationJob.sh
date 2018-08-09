#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source P.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
source Jellyfish.sh
FSQ=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.fq" } else { print "/share/hormozdiarilab/Codes/NebulousSerendipity/data/" GEN ".Fq" } }')
echo $FSQ
echo "$@"
python -m kmer.depth --job UniqueKmersDepthOfCoverageEstimationJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --jellyfish $JLY $RJF --fastq $FSQ --ksize $KSZ $SIM "$@"
