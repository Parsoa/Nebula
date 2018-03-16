#!/bin/bash

#SBATCH --partition=gc
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --job-name="HG512/EstimateDepthOfCoverage"
#SBATCH --ntasks=24
#SBATCH --time=3-12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/HG512.coverage.out
#SBATCH --mail-type=ALL

cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv/bin/activate
cd src/python
python -m kmer.prune --threads 24 --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/hg38.exons.filtered.bed --ref hg38 --fastq /share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio_Phase4/CHS_Trio/HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.Test.Fq
