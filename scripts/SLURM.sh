#!/bin/bash
export P=$(pwd)
source BED.sh
echo PWD $P
JOB=$(echo $P | awk -F/ '{ print $NF }')
echo JOB $JOB
echo MEM $1
echo TME $2
touch "$P"/submit.sh
chmod +x $P/submit.sh
# SLURM Headers
echo "#!/bin/bash" > $P/submit.sh
echo "#SBATCH --partition=production" >> $P/submit.sh
echo "#SBATCH --nodes=1" >> $P/submit.sh
echo "#SBATCH --mem=$1" >> $P/submit.sh
echo "#SBATCH --job-name=$BED.$GEN.$JOB" >> $P/submit.sh
echo "#SBATCH --ntasks=48" >> $P/submit.sh
echo "#SBATCH --time=$2" >> $P/submit.sh
echo "#SBATCH --mail-user=pkhorsand@ucdavis.edu" >> $P/submit.sh
echo "#SBATCH --output=$P/slurm.out" >> $P/submit.sh
echo "#SBATCH --mail-type=ALL" >> $P/submit.sh
# Code
echo 'export PYTHONPATH=$PYTHONPATH:/home/pkhorsand/local/cplex/lib/python' >> $P/submit.sh
echo 'source PWD.sh' >> $P/submit.sh
echo 'source SIM.sh' >> $P/submit.sh
echo 'cd /share/hormozdiarilab/Codes/NebulousSerendipity' >> $P/submit.sh
echo 'source venv2/bin/activate' >> $P/submit.sh
echo 'cd src/python' >> $P/submit.sh
echo 'source BED.sh' >> $P/submit.sh
echo 'source FASTQ.sh' >> $P/submit.sh
echo 'source JELLYFISH.sh' >> $P/submit.sh
echo 'echo "$@"' >> $P/submit.sh
shift 2
echo 'python -m kmer.main --job $JOB --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/$BED --threads 48 --reference $REF --jellyfish $JLY $RJF --fastq $FSQ --genome $GEN --readlength 100 $SIM --seed 165784623 --heterozygous $DESCP' $@ >> $P/submit.sh
