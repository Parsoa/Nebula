#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --job-name="prune"
#SBATCH --ntasks=48
#SBATCH --time=0-01:30:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python

python -m kmer.prune --reference hg19 --bed CHM1_Lumpy.Del.100bp.DEL.bed --fastq /share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq
