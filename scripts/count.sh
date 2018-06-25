#!/bin/bash
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
python -m kmer.genotyping --job LocalUniqueKmersCountingJob --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/HG00512.hg38.DEL.bed --threads 48 --fastq /share/hormozdiarilab/Codes/NebulousSerendipity/HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.Test.ReSorted.Fq
