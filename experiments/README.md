# Introduction

This README files provides instructions on how to reproduce the results presented on the journal version of Nebula, submitted to Oxford Nucleic Acid Research in June 2020.

# Data

We used the SV calls made by <cite>[Chaisson et al][1]</cite> on 1KG samples [HG00514](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/HG00514/), [HG00733](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/PUR/HG00733/) and [NA19240](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/YRI/NA19240/high_cov_alignment/) as the basis of this study. The deletion and insertion calls can be retrievd from [this link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180627_PanTechnologyIntegrationSet/). The inversion calls can be retrived from dbVar's ftp website [here](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/).

## Input Preparation

As Nebula genotypes SVs purely based on coordinates, we have done slight modifications to the calls above to make genotyping results more consistent.

* For overlapping deletions, we only keep the one with the smallest `BEGIN` position and discard the rest.

* Events that are too close to each other will yield more or less the same set of kmers and can be genotyped in place of one another. We merge all insertions reported withing 50 bases of each other into a single insertion across all three samples for consistency.

* We have removed events from repeat regions, i.e those with `IS_TRF=True`.

Download the above VCF files, then use the Python scripts to perform the filtering:

```
cd /nebula/clone/directory
virtualenv -p python3 venv
source venv3/bin/activate
pip install -r requirements.txt
# download VCF files and store somewhere
cd src/python
python -m nebula unify --vcf HG00514.merged_nonredundant.vcf HG00733.merged_nonredundant.vcf NA19240.merged_nonredundant.vcf --workdir /directory/to/save/modified/VCF/files
```

Use the `parse_vcf.py` script to convert the VCF files into BED files for input to Nebula:

```
# all variants
./parse_vcf.py HG00514.merged_nonredundant.unified.vcf > HG00514.unified.bed
./parse_vcf.py HG00733.merged_nonredundant.unified.vcf > HG00733.unified.bed
./parse_vcf.py NA19240.merged_nonredundant.unified.vcf > NA19240.unified.bed
```

# Running

We now provide instruction on how to reproduce the results on NA19240 presented in the study.

First, download the set of kmers for estimating sample coverage and GC correction from the links below:

** add links

Now we can extract kmers from HG00514 and HG00733:

```
./scripts/preprocess.sh
```

This extracts the kmers and stores them in the directory `kmers`. This step should take ~15 minutes. Now we need to genotype both samples with the exrtacted kmers and filter the potentially misleaing kmers:

```
./scripts/genotype_HG00514.sh
./scripts/genotype_HG00733.sh
```

This should take around 1 hour for each sample. Samples can be genotyped simultaneously.

Finally, the kmers are to be filtered and exported for genotyping:

```
./scripts/mix.sh
```

This stores the kmers in `Mix`.

Now we can genotype NA19240:

```
./scripts/genotype_NA19240.sh
```

The output file `genotypes.bed` under the directory `NA19240` contains the genotype predictions.

You can adjust the paths in the scripts according to where the various files are located. All paths must be absolute.

