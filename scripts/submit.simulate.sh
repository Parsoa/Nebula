#!/bin/bash

#SBATCH --partition=intel
#SBATCH --nodes=1
#SBATCH --mem=160G
#SBATCH --job-name="Simulation"
#SBATCH --ntasks=24
#SBATCH --time=12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/simulation/HG00513.hg38.chr1.DEL.bed/31/Simulation/simulation.out
#SBATCH --mail-type=ALL
cd /share/hormozdiarilab/Codes/NebulousSerendipity/simulation/HG00513.hg38.chr1.DEL.bed/31/Simulation
/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/simulate.sh
