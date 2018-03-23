#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --job-name="CHM1/MostLikelyBreakPoints"
#SBATCH --ntasks=24
#SBATCH --time=12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/CHM1.MostLikelyBreakPoints.out
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python
python -m kmer.break_point --coverage 50 --std 20 --threads 24
