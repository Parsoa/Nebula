#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --job-name="prune"
#SBATCH --ntasks=24
#SBATCH --time=6:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python
./run.sh
#python -m kmer.prune --fastq /share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq --threads 48
#python -m kmer.prune --reference hg38 --bed /share/hormozdiarilab/Data/ReferenceGenomes/Hg38/wgEncodeGencodeExonSupportV27.txt --fastq /share/hormozdiarilab/Data/Genomes/Illumina/CHMs/CHM1_hg38/CHM1.samtoolsversion.fq --threads 
python -m kmer.break_point --coverage 30 --std 20 --threads 244
