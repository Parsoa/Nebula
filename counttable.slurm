#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=70G
#SBATCH --job-name="prune"
#SBATCH --ntasks=8
#SBATCH --time=4-12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python

#python -m kmer.prune --fastq /share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq --threads 48
python -m kmer.count_server --threads 4 --counttable /share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio/HG00512.cor.fq
