#!/bin/bash
export P=$(pwd)
echo PWD $P
JOB=$(echo $P | awk -F/ '{ print $NF }')
echo JOB $JOB
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source P.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
source FASTQ.sh
source Jellyfish.sh
echo "$@"
python -m kmer.main --job $JOB --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --reference $REF --jellyfish $JLY $RJF --fastq $FSQ --readlength 100 --insertsize 500 --gap 5 --ksize $KSZ $SIM $DESCP "$@"
#/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/$JOB.sh $@
