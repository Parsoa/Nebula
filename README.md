# Nebula

Nebula is an ultra-efficient mapping-free structural variation genotyper based on kmer counting.

# Operating Principles

Nebula is a mapping-free approach for accurate and efficient genotyping of SVs. Nebula is a two-stage approach and consists of a **kmer extraction phase** and a **genotyping phase**. Given as input a set of SV coordinates (BED/VCF), the  reference assembly (FASTA), and a set of mapped samples on which the genotype of these SVs is already known, Nebula extracts a collection of kmers that represent the input SVs (kmer extraction phase) and these extracted kmers will then be used to genotype the same set of SVs on any new WGS sample(s) without the need to map the reads to the reference genome (genotyping phase). This is done by counting the kmers in the WGS reads of the new sample(s) and predicting genotypes using a likelihood model. 

![Nebula's pipeline](assets/Pipeline.png)

# Installation

Nebulo requires htslib to be installed for parsing BAM and SAM files but is otherwise self-contained.

Download the latest executable from releases section and put it somewhere inside your PATH.

# Usage

## kmer Extraction 

For kmer extraction, a set of BAM/SAM/CRAM files should be passed to Nebula along with one BED file with the genotypes of the target SVs on that sample for each BAM file.
The BED file can minimally contain the event coordinates and types: 
```
#CHROM BEGIN   END SVTYPE
```

The `SVTYPE` should be one of `DEL`, `INS` or `INV`. Optionally, the length of the SV and the actual inserted sequence can be passed for insertions (`SVLEN` and `SEQ`). Any additional field is ignored. 

Nebula expects a certain directory structure for outputs of different stages. 

```
nebula preprocess --bed /path/to/genotypes_1.bed /path/to/genotypes_2.bed --bam /path/to/bam_file_1.bed /path/to/bam_file_2.bed --wokdir output/kmers --reference /path/to/reference/FASTA/file --thread <number of threads to use>
```

This will output a number of JSON files including the kmers in the direcotry set by `workdir`.

Next, the input samples should be genotyped with these kmers. Create a directory for this stage's output, e.g `output`. The genotyping output for each of the samples must be stored in subdirectory inside `outout` with the same name as the sample. A sample's name is just whatever identificationn you use for that sample, but has to consistent through the pipeline:

```
nebula genotype --bed /path_to_genotypes_1.bed --bam /path/to/bam_file_1.bed --workdir output/sample_1 --kmers /output/kmers --depth_kmers depth_kmers.json --gc_kmers gc_kmers.json
nebula genotype --bed /path_to_genotypes_2.bed --bam /path/to/bam_file_2.bed --workdir output/sample_2 --kmers /output/kmers --depth_kmers depth_kmers.json --gc_kmers gc_kmers.json
```

Merge the remaining kmers after filtering. Note that this stage will determine the output directory for each sample based on the workdir and the name of each sample: 

```
nebula mix --bed /path_to_genotypes_1.bed//path_to_genotypes_2.bed --samples sample_1,sample_2 --workdir ./output
```

The output kmers are stored in afolder named `Mix` inside workdir.

# Genotyping

For genotyping unmapped sample with the extracted kmers from an earlier kmer-extraction run:

```
nebula.sh genotype --kmers /path/to/Mix/directory --bam/--fastq /path/to/sample --workdir <output directory>
```

Nebula will output a BED file named `genotypes.bed` in the specified working directory. The file will include the original fields in the input BED files along with the field `GENOTYPE` (one of 0/0, 1/0 or 1/1). Note that a BED file does not need to passed to the genotyper; the variants are implicit in the kmers. There are no requirements on the output directory.

# Benchmarking and Performance

Nebula is designed to be simple, fast and memory efficient so it can be run on any reasonable personal hardware. Using a single processor core, Nebula can count kmers at a rate of 400,000 reads per second from a FASTQ file. A 30x human sample can be process in less than 80 minutes on a single core.

# Citation

The pre-print is currently available on BioRxiv but it's very outdated:

The current jounral version has been completely rewritten and is now under review.

