#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH --job-name="CHM1/MostLikelyBreakPoints"
#SBATCH --ntasks=48
#SBATCH --time=12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/CHM1.most_likely.out
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python
python -m kmer.break_point --coverage 30 --std 9 --threads 24
