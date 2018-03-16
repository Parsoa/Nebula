#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --job-name="HG512/EstimateDepthOfCoverage"
#SBATCH --ntasks=24
#SBATCH --time=12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/CHM1.BreakPointJob.out
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python
python -m kmer.break_point --threads 24
