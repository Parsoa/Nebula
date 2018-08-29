#!/bin/bash
#SBATCH --partition=intel
#SBATCH --nodes=1
#SBATCH --mem=90G
#SBATCH --job-name="NebulousSerendipity"
#SBATCH --ntasks=48
#SBATCH --time=12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/current.out
#SBATCH --mail-type=ALL
P=$(pwd)
echo $P
CMD = $(echo $P | awk -F/ '{ print $NF }')
srun --mem=12000 --mincpus=48 --time=6000
