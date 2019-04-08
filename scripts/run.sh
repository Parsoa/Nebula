#!/bin/bash
export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python
source PWD.sh
source SIM.sh
source CMD.sh
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
source BED.sh
source BAM.sh
source FASTQ.sh
source JELLYFISH.sh
echo "$@"
python -m kmer.main $CMD --job $JOB --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --reference $REF --jellyfish $RJF --fastq $FSQ --bam $BAM --genome $GEN --readlength 100 $SIM --seed 165784623 $DESCP --workdir /share/hormozdiarilab/Codes/NebulousSerendipity/output "$@"
