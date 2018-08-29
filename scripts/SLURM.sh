#!/bin/bash
export P=$(pwd)
echo PWD $P
JOB=$(echo $P | awk -F/ '{ print $NF }')
echo JOB $JOB
echo MEM $1
echo TME $2
touch "$P"/submit.sh
chmod +x "$P"/submit.sh
echo "#!/bin/bash" > $P/submit.sh
echo "#SBATCH --partition=production" >> $P/submit.sh
echo "#SBATCH --nodes=1" >> $P/submit.sh
echo "#SBATCH --mem=$1" >> $P/submit.sh
echo "#SBATCH --job-name=$JOB" >> $P/submit.sh
echo "#SBATCH --ntasks=48" >> $P/submit.sh
echo "#SBATCH --time=$2" >> $P/submit.sh
echo "#SBATCH --mail-user=pkhorsand@ucdavis.edu" >> $P/submit.sh
echo "#SBATCH --output=$P/slurm.out" >> $P/submit.sh
echo "#SBATCH --mail-type=ALL" >> $P/submit.sh
cat /share/hormozdiarilab/Codes/NebulousSerendipity/scripts/$JOB.sh >> $P/submit.sh
