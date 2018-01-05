#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --job-name="prune"
#SBATCH --ntasks=2
#SBATCH --time=1-12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python

#python -m kmer.prune --fastq /share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq --threads 48
python -m kmer.prune --reference hg38 --bed /share/hormozdiarilab/Data/ReferenceGenomes/Hg38/wgEncodeGencodeExonSupportV27.txt --fastq /share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq --threads 1
