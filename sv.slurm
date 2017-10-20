#!/bin/bash

#SBATCH --partition=gc # partition to submit to
#SBATCH --job-name=“nebula_sv”
#SBATCH --nodes=1 # single node
#SBATCH --ntasks=4 #
#SBATCH --mem=80G # memory pool all cores, default is 2GB per cpu
#SBATCH --time=08:00:00 # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/output/chr6_130848185_130852294.sv # STDOUT
#SBATCH --error=/share/hormozdiarilab/Codes/NebulousSerendipity/output/chr6_130848185_130852294.sv.err # STDOUT
#SBATCH --mail-user=pkhorsand@ucdavis.edu

/share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd ./src/python
python -m kmer.count_server
python -m kmer.sv