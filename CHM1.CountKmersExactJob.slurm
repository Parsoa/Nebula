#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --job-name="CHM1/CountKmersExactJob"
#SBATCH --ntasks=1
#SBATCH --time=1-12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/CHM1.CountKmersExactJob.out
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/cpp
#g++ -std=c++11 merge_kmer_counts.cpp -static-libstdc++
./a.out
