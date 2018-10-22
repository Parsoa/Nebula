#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source PWD.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
source FASTQ.sh
source JELLYFISH.sh
echo "$@"
python -m kmer.main --job $JOB --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --reference $REF --jellyfish $JLY $RJF --fastq $FSQ --genome $GEN --readlength 100 --insertsize 500 --gap 5 --ksize $KSZ $SIM --seed 165784623 --heterozygous $DESCP "$@"
