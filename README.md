# Nebula

Nebula is a ultra-efficient tool for genotyping structural variations without the use of read mapping.

# Installation

Nebula depends on https://Github.com/coin-or/Clp for solving the linear programming model. Several other common bioinformatis libraries are also required for the 
preprocessing stage. There are scripts provided for installing all dependencies.

Clone this repository and run:

`
virtualenv venv
pip install -r requirements.txt
./htslib.sh
./jellyfish.sh
./coin.sh
./counter.sh
`

Then add `/src/python/nebula` to PYTHONPATH or alternatively add the `scripts` directory to PATH.

# Basic Usage

Nebula runs in two steps:

1. Preprocessing
2. Genotyping

## Preprocessing

This stage requires an input BED file with a set if structural variation coordinates. Currently deletions, insertions
and inversions are supported.

The preprocessor can be run as below:

`python -m nebula --bed <path to bed file> --ref <path to a reference assembly> --workdir <working directory> --threads <number of threads to use> --bam <path to a bam file with structural variations>`

The BED file should include the structural variation tracks with the standard header:

`CHROM BEGIN END SVTYPE`

Where `SVTYPE` is one of `DEL`, `INS`, `INV`, `ALU` or `MEI`. If not provided, a default of `DEL` will be assumed.

Events inserting new sequences such as insertions require a `SEQ` field with the inserted sequence.

An optional `ID` field is supported. This will be the id used internally and in output files for the event. If non provided, a feault of `CHROM_BEGIN_END_SVTPYE` will be used instead.

Ideally, this BED file includes events previously genotyped on a specific sample by (most likely) a mapping-based genotyper. If that's the case, it is recommended to
provide the mapped BAM file to Nebula. This enables the use of **junction kmers** that can hugely improve genotyping accuracy.

The preprocessor will create a number of directories inside `workdir`. The kmers signatures for the given events will be available in
`workdir/MixKmersJob/kmers.json`. These JSON file can be fed into the genotyper in the next steps.

# Genotyping

Once kmers are extracted, any sample can be genotyped. The BED file and its corrsponding set of kmers should be propvided to the genotyper as arguments.
Multiple sets of kmers from independent BED files can be passed. Nebula will merge kmers, remove duplicates and report genotypes for the union of all events provided. 

Nebula only requires the raw set of reads in FASTQ format genotyping. Reads stored in BAM format are also supported.

To genotype:

`python -m nebula genotype --bed <set of bed files> --kmers <kmers.json for each preprocessed bed file> --bam/--fastq <path to sequencing reads> --threads <default 1> --workder <output directory>`

Most of the genotyper's runtime is spent counting kmers in the set of reads. This can be parallelized if a FASTQ files is used but is not supported for BAM files. If enough memory bandwidth is available, Nebula will
grow linearly faster with additional number of threads. The LP solve stage can also benefit from more threads but is fast enough on its own.

Nebula will output a BED file named `genotypes.bed` in its working directory. The genotyping and LP values plus certainty scores will be included as new columns
`GENOTYPE`, `LP_VALUE`, `SCORE`.

# Clustering

Nebula can optionally cluster genotypes on multiple samples to potentially improve results. This is done via

`python -m nebula cluster --bed <the output bed files for different samples> --workdir`

This clusters based on the LP value and doesn't look at the reads or kmers counts. Only the intersection of all provided BED files will be clustered, so this mode is 
meant to be used on genotyping output of the same set of kmers.

# System Requirements

Nebula is designed to be very fast and memory efficient so it can be virtually run on any reasonable personal hardware. A FASTQ files with 2 billion reads can be processed in about 2 hours on a single Core i7 core.

Nebula will require around 1GB of memory for genotyping 10,000 SVs. This will grow linearly with number of threads, so it is advised to generally run with a single core on low-end hardware. 
On clusters, it may be more beneficial to genotype many samples concurrently rather than using many CPU threads to genotype a single sample.
