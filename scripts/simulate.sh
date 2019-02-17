#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source PWD.sh
source SIM.sh
source KSIZE.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
source BAM.sh
source FASTQ.sh
source JELLYFISH.sh
echo "$@"
python -m kmer.main --job $JOB --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --reference $REF --jellyfish $JLY $RJF --fastq $FSQ --bam $BAM --genome $GEN --readlength 100 --gap 5 --ksize $KSZ $SIM --seed 165784623 $DESCP "$@"
cd $P
cat chr*.1.fq > test.1.fq &
cat chr*.2.fq > test.2.fq &
