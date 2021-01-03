# Introduction

This README files provides instructions on how to reproduce the results presented on the journal version of Nebula, submitted to Oxford Nucleic Acid Research in June 2020.

# Data

We used the SV calls made by [Chaisson et al](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180627_PanTechnologyIntegrationSet/) on 1KG samples [HG00514](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/HG00514/), [HG00733](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/PUR/HG00733/) and [NA19240](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/YRI/NA19240/high_cov_alignment/) as the basis of this study. The deletion and insertion calls can be retrievd from [this link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180627_PanTechnologyIntegrationSet/). The inversion calls can be retrived from dbVar's ftp website [here](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/).

For consistency, we have merged events that overladp or are less than 10bp aprt (only for insertions) into single events across all three samples. This allows for a more streamlined comparison of genotypes. SVs that are form repeat regions (`IS_TRF=TR`) and SVs shorter than 10bp have been filtered. Insertions with an absent `SEQ` field have also been filtered. The original SV calls and the modified calls can be found in the `variants` folder.

# Running

The set of kmers for estimating sample coverage and GC coverage can be found in the `kmers` folder. The genotyping kmers have to be generated locally as they are too large (1GB on disk) to include in Github.

You will need to download the BAM files for HG00514 and HG00733. For NA19240 only a FASTQ file is needed, although Nebula can also genotype a BAM file.

To genotype the SVs from non-repeat regions on HG00514 and HG00733 on NA19240, follow the steps below. All commands are run from this directory.

1. Extract genotyping kmers from HG00514 and HG00733:

```
nebula preprocess --bed $PWD/variants/HG00514.unified.bed,$PWD/variants/HG00733.unified.bed --reference GRC38.fasta --bam /path/to/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam,/path/to/HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.bam --workdir $PWD/kmers --threads 4
```

This extracts the kmers and stores them in the directory `kmers`. This step should take ~20 minutes.


2. Genotype both samples with the exrtacted kmers and filter the potentially misleaing kmers:

```
nebula genotype --bed $PWD/variants/HG00514.unified.bed --gc_kmers $PWD/kmers/gc_kmers.json --depth_kmers $PWD/kmers/reference_kmers.json --bam /path/to/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam --kmers $PWD/kmers --workdir $PWD/HG00514 --threads 8 --select --unique
nebula genotype --bed $PWD/variants/HG00733.unified.bed --gc_kmers $PWD/kmers/gc_kmers.json --depth_kmers $PWD/kmers/reference_kmers.json --bam /path/to/HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.bam --kmers $PWD/kmers --workdir $PWD/HG00733 --threads 8 --select --unique
```

This should take around ~70 minutes for each sample. Samples can be genotyped simultaneously.

3. Merge the filtered kmers to get the combined set of kmers for HG00514 and HG00733:

```
nebula mix --bed $PWD/variants/HG00514.unified.bed,$PWD/variants/HG00733.unified.bed --samples HG00514,HG00733 --workdir $PWD --threads 8
```

4. Now we can genotype NA19240:

```
nebula genotype --bed $PWD/variants/HG00514.unified.bed,$PWD/variants/HG00733.unified.bed --gc_kmers $PWD/kmers/gc_kmers.json --depth_kmers $PWD/kmers/reference_kmers.json --fastq /path/toNA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.bam --kmers $PWD/Mix --workdir $PWD/NA19240 --threads 8
```

The output file `genotypes.bed` under the directory `NA19240` contains the genotype predictions.

Nebula requires absolute paths for all arguments. If you are inside a different directory, change the paths as needed.

