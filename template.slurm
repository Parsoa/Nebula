#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=90G
#SBATCH --job-name="PairedEndCounting"
#SBATCH --ntasks=48
#SBATCH --time=4-12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/HG512.coverage.out
#SBATCH --mail-type=ALL
cd 
./run.sh
1
2
