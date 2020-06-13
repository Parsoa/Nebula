# Introduction

This README files provides instructions on how to reproduce the results presented on the journal version of Nebula, submitted to Oxford Nucleic Acid Research in June 2020.

# Data

We used the SV calls made by <cite>[Chaisson et al][1]</cite> on 1KG samples [HG00514](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/HG00514/), [HG00733](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/PUR/HG00733/) and [NA19240](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/YRI/NA19240/high_cov_alignment/) as the basis of this study. The deletion and insertion calls can be retrievd from [this link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180627_PanTechnologyIntegrationSet/). The inversion calls can be retrived from dbVar's ftp website [here](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/).

## Deletions and Insertions

As Nebula purely cares about the coordinates of the SV breakpoints, we have done slight modifications to the calls above to make genotyping results more consistent.

For overlapping deletions, we only keep the one with the smallest BEGIN position and discard the others.

For insertions, we are assuming that the `END` field in the `INFO` section is set to `BEGIN + 1`. We have modified calls that don't adhere to this assumption. Sometimes, the same insertion is reported on slightly different coordiates between the three samples. For insertions less than 100bp apart, we are unifying all into the one with the smallest `BEGIN` position. To run the unification step, download the VCF files from the link above and run Nebula as below:

```
nebula.sh unify --vcf <path to VCF files> --workdir <directory to store unified BED files>
```

This produces a set of unified VCF files for each sample. For instance `HG00514.merged_nonredundant.unified.all.vcf` includes all events on HG00514 while `HG00514.merged_nonredundant.unified.repeat.vcf` only includes events in repeat regions and `HG00514.merged_nonredundant.unified.vcf` includes only non-repeat events. Repeat events are identified by the presence of `IS_TRF=True` in the HGSV VCF files. We only used the non-repeat events for performance evaluation. It is preferable to pass BED files to Nebula as input. The VCF files can be converted to BED format using the `parse_vcf.py` script (note the `-` in arguments):

```
./scripts/parse_vcf.py HG00514.merged_nonredundant.vcf - > HG00514.merged_nonredundant.bed
```

The other tools used in this study may require additional VCF fields not present in the data above, we have provided trivial values for those fields. There is a script that transforms the input VCF file for the required format for each tool. To convert a VCF file into Paragraph's required format run:

```
paragraphify_vcf.py HG00514.merged_nonredundant.unified.vcf
# outputs a file HG00514.merged_nonredundant.unified.paragraph.vcf
```

The VCF files for HG00514 and HG00733 need to be merged as most tools only accept one VCF files for input (Nebula can take multiple BED files as input and will remove duplicates and merge overlapping events internally). The script `merge.sh` can be used for this purpose. It requires the presence of `HG00514.merged_nonredundant.unified.paragraph.vcf` and `HG00733.merged_nonredundant.unified.paragraph.vcf`. The scripts outputs two files `HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.DEL.vcf` with deletions and `HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.INS.vcf` with only insertions (the separation is necessary as Delly and SVTyper can't genotype insertions). These files can be used for all tools except BayesTyper.

## Inversions

For inversions we have a single VCF file that includes genotypes for all 9 samples in CHS, PUR and YRI trios (`HGSV_ILL.Integration.revise.20180208.vcf`). This file includes deletions and insertions as well as other SVs. To only get the inversions run:

```
grep -E "#|INV" HGSV_ILL.Integration.revise.20180208.vcf >  HGSV_ILL.Integration.revise.20180208.INV.vcf
```

Run `parse_inv_vcf.py` on this file to get a set of BED files for each of the three samples:

```
./parse_inv_vcf.py HGSV_ILL.Integration.revise.20180208.INV.vcf
```

For the input to other tools, run the `paragraphify_vcf.py` script:

```
./paragraphify_vcf.py HGSV_ILL.Integration.revise.20180208.INV.vcf
```

The output can be used by SVtyper and Paragraph.

# Running

After the unification step, we have the same of events across all three samples. This makes the comparison results more consistent and reproducable.

## Nebula

Nebula first needs to preprocess HG00514 and HG00733 for kmers. The kmers can be downloaded from [here](). Alternatively, th kmers can be reproduced locally by running the preprocessing stage. Before running the preprocessing stage one needs to extract kmers for GC content estimation for Nebula. This stage selects kmers from different regions of the reference genome with different GC content levels. These kmers will be counted in the genotyping sample to provide better estimates of coverage across the genome. This need only be done once and the kmers are only dependent on the reference.

```
# Extract GC content esimtation kmers
nebula.sh gc --workdir output/preprocessing/GC_Kmers --reference GRC38.fasta --jellyfish hg38_mer_counts_32k.jf

# For HG00514
nebula.sh preprocess --bed HG00514.unified.bed --reference GRC38.fasta --bam HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam --workdir output/preprocessing/HG00514 --jellyfish hg38_mer_counts_32k.jf --gckmers gc_kmers.json
nebula.sh genotype --bed HG00514.unified.bed --reference GRC38.fasta --bam HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam --workdir output/preprocessing/HG00514 --kmers output/preprocessing/HG00514/MixKmersJob/kmers.json --select
# HG00733
nebula.sh preprocess --bed HG00733.unified.bed --reference GRC38.fasta --bam HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.bam --workdir output/preprocessing/HG00733 --jellyfish hg38_mer_counts_32k.jf --gckmers gc_kmers.json
nebula.sh genotype --bed HG00733.unified.bed --reference GRC38.fasta --bam HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.bam --workdir output/preprocessing/HG00733 --kmers output/preprocessing/HG00733/MixKmersJob/kmers.json --select
```

This should take around 2 hours for each sample. 

We ran the inversion experiments later and separately, so we repeated the above steps for inversions. Alternatively, once can pass BED files for inversions along with the one for deletions and insertions and get the same results. This avoids the overhead of running the preprocessing for the same sample twice.

Once completed, the third sample can be genotyped. First (optionally) convert the BAM file for NA19240 into FASTQ using `bamtofastq`, then run:

```
nebula.sh genotype --bed HG00514.unified.bed HG00733.unified.bed --bam NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.fastq --workdir output/genotyping/NA19240/HG00514_HG00733 --kmers output/preprocessing/HG00514/ExportGenotypingKmersJob/kmers.json output/preprocessing/HG00733/ExportGenotypingKmersJob/kmers.json
```

One can also pass BAM file to Nebula directly for genotyping, there is no advantage in this as the approach is mapping-free, and counting kmers in BAM files will be slightly slower than FASTQ files.

## BayesTyper

BayesTyper has a complex pipeline for selecting variants to genotype, however we are only inputting the unified 1KG variants here and the pipeline can be skipped. BayesTyper needs to preprocess the input VCF file and requires KMC tables for genotyping. We use the same merged VCF file as before for deletion and insertions.

First download the data bundle from [BayesTyper's Github repo](https://github.com/bioinformatics-centre/BayesTyper) for GRCh38 and store in `output/genotyping/NA19240/HG00514_HG00733/Verification/BayesTyper`

Navigate to this directory and create a file named `samples.tsv` with the contents below:

```
NA19240	F	<path to NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.bam>
```

Create the KMC index for the sample:

```
kmc -k55 -ci1 NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.Fq <output_prefix> ./KMC
```

Next create Bloom Filters for BayesTyper:

```
bayesTyperTools makeBloom -k <output_prefix> -p 16
```

Next, the VCF file alleles have to converted to BayesTyper's required format:

```
bayesTyperTools convertAllele -v HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.vcf -g BayesTyper/GRCh38.fa --keep-imprecise -o HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.bayestyper.vcf
```

Now BayesTyper needs to cluster the variants:

```
bayesTyper cluster -v HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.bayestyper.vcf -s samples.tsv -g BayesTyper/GRCh38_canon.fa -d BayesTyper/GRCh38_decoy.fa -p 16
```

Now BayesTyper can genotype the calls:

```
bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data -s samples.tsv -g BayesTyper/GRCh38_canon.fa -d BayesTyper/GRCh38_decoy.fa -o bayestyper_unit_1/bayestyper -z -p 4 --disable-observed-kmers --noise-genotyping
```

This outputs a file named `bayestyper.vcf` in the directory `bayestyper_unit_1` under the current working directory.

## Paragraph

Create the directory `Paragraph` under `output/genotyping/NA19240/HG00514_HG00733/Verification`.

Download and compile Paragraph as per the instruction in their Github repo. Create the file `samples.txt`:

```
#id,path,depth,read length
NA19240,NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.bam,40,100
```

Run with default options:

```
<path to multigrmpy.py> -r GRC38.fasta -i HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.DEL.vcf -m samples.txt -o ./DEL
```

This will genotype the deletions and output a file named `genotypes.vcf` under the `DEL` directory.

## SVTyper

Create a directory for SVTyper similar to other tools. Install the package and run with default options:

```
svtyper -i HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.DEL.vcf -B NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.bam > DEL/genotypes.vcf
```

## Delly

Create a directory for Delly's output and run with default options. Note that Delly requires a reference with deocy sequences and a file with regions to be exluded (available on their Github repo):

```
delly call -g GRCh38_decoys.fa -v HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.DEL.vcf -x human.hg38.excl -o ./DEL/genotypes.bcf NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.bam
```

# Validation

Nebula outputs a file named `merge.bed` unde the directory `output/genotyping/NA19240/HG00514_HG00733/CgcIntegerProgrammingJob` which includes the predicted genotypes. To compare Nebula's results against actual NA19240 predictions  run:

```
tabulate.sh merge.bed
verify.sh NA19240.unified.bed
```

This creates several files:

```
00_as_00.bed
00_as_10.bed
00_as_11.bed
10_as_00.bed
10_as_10.bed
10_as_11.bed
11_as_00.bed
11_as_10.bed
11_as_11.bed
```

With 3 genotypes of 0/0, 0/1 and 1/1, there are 9 possibilities for each prediction (e.g a 0/0 event being predicted as 1/0, etc). Each file above contains the events falling in one of these 9 categories.

The script `intersect.sh` works similar to `bedtools intersect` however is specifically suited to our validation process.

To compare Nebula's results with those of Paragraph, convert Paragraph's output to BED format using the `parse_paragraph_vcf.py` script. Then run `verify_paragraph.sh`. This outputs genotype statistics for Paragraph and Nebula. Similarly, repeat the procedure for BayesTyper, Delly and SVtyper to get their results to get the number for each tool using the corresponding script.

For the comparison, we have only considered events that one of Delly, SVtyper or Paragraph can correctly predict as present (exact genotype of 1/0 or 1/1 does not matter) on HG00514 and HG00733. This requires us to run all tools on HG00514 and HG00733 first and use th set of events predicted as present as the gold standard. The scripts used above assume these events are stored in files `consensus.DEL.bed`, `consensus.INS.bed` and `consensus.INV.bed`.

# Citations

[1]: https://www.nature.com/articles/s41467-018-08148-z

