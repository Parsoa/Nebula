# Nebula

Nebula is an ultra-efficient mapping-free structural variation genotyper based on kmer counting.

# Operating Principles

Nebula works in two stages. During the preprocessing phase, Nebula takes as input a mapped WGS sample and a set of SVs present on this sample and extracts kmers whose frequency in the sequencing reads is dependent on the presence of said SVs. These kmers can then be used to genotype the same set of structural variations on any other sample. The input to this stage is usually the VCF output of a discovery tool such a Delly or Lumpy.

During the genotyping phase, Nebula uses the kmers from a previous preprocessing stage to genotype a new sample. Multiple sets of kmers can be passed at the same time and the union of events represented by them will be genotyped. Kmers are counted in the input FASTQ files and the counts are plugged into a linear programming model that outputs the final genotypes after rounding. This stage is mapping-free, as a result large libraries of newly sequenced unmapped samples can be efficiently genotyped with Nebula.

![Nebula's pipeline](assets/Pipeline.png)

# Installation

Nebula depends on [Clp](https://Github.com/coin-or/Clp) for solving the linear programming model and on [Jellyfish](https://github.com/gmarcais/Jellyfish) for the preprocessing stage. There are scripts provided for installing all dependencies. Nebula requires Python 3.5 and higher.

Clone this repository and run:

```
virtualenv venv
pip install -r requirements.txt
./htslib.sh # install htslib
./jellyfish.sh # install Jellyfish
./coin.sh # install COIN-OR
./counter.sh # compile sources
```

Finally copy the `nebula.sh` script to a location in `PATH`.

# Usage

## Preprocessing

Multiple BED/VCF files can be passed as input to the preprocessing stage together however only one sample can be preprocessed at a time. All input files should contain events on the same sentence and from the same reference version.

The preprocessing stage requires a reference assembly in FASTA format and the Jellyfish index for that assembly with kmer length 32. The index need only be created once and can be reused on all future runs:

```
jellyfish index
```

Nebula requires to extract a set of kmers from regions with different GC content levels across the genome. This is done once and the resulting kmers can be reused for any future runs:

```
nebula.sh gc --workdir <where to store kmers> --reference GRCh38.fa --jellyfish hg38_mer_counts_32k.jf
```

Alternatively, a pre-computed Jellyfish index for GRCh38 and set of GC kmers can be downloaded from this repository. In that case, the `--jellyfish` option can be ommitted below and there is no need to install Jellyfish.

Samples usage for the preprocessor is as follows:

```
nebula.sh preprocess --bed <path to bed file> --reference <path to a reference assembly> --workdir <working directory> --bam <path to a bam file with structural variations> --jellyfish <path to Jellyfish index for the reference> --gckmers <path to gc_kmers.json>
```

The BED file can minimally contain the event coordinates and types: 

`#CHROM BEGIN   END SVTYPE`

The `SVTYPE` should be one of `DEL`, `INS` or `INV`.

Optionally, the length of the SV and the actual inserted sequence can be passed for insertions (`SVLEN` and `SEQ`). The fiels will be ignored for deletions and inversions. An optional `ID` field is supported. This will be the id used internally and in output files for the event. If not provided, a default value of `SVTYPE@CHROM_BEGIN_END_SVTPYE` will be used instead. Any additional field is ignored. 

The preprocessor can optionally use the extracted kmers to genotype the given sample and will only retain the events that could be genotyped as `1/0` or `1/1` (with `--select`). This may reduce the number of events that can be genotyped on other samples (recall) but should increase the precision of those genotype calls. Kmers will be output under `ExportGenotypingKmersJob/kmers.json`. A single set of kmers will be exported no matter the number of input BED files.

# Genotyping

The genotyper takes as input one or more sets of kmers from previous preprocessing runs along with their corresponding BED tracks and a FASTQ file to be genotyped. Multiple FASTQ files can be passed but they will be considered to be parts of the same sample. Alternatively, a single BAM files can be passed, however only the raw read sequences are used.

```
nebula.sh genotype --bed <set of bed files> --kmers <kmers.json for each preprocessed bed file> --bam/--fastq <path to sequencing reads> --workdir <output directory>
```

Nebula will output a BED file named `merge.bed` in its working directory under `CgcIntegerProgrammingJob`. The file will include the original fields in the input BED files along with the fields `GENOTYPE` (one of 0/0, 1/0 or 1/1), `LP_VALUE` which is the real-valued genotype from the linear program and `CONFIDENCE` which can be either HIGH or LOW, representing Nebula's confidence in the correctness of the genotype.

# Benchmarking and Performance

Nebula is designed to be simple, fast and memory efficient so it can be run on any reasonable personal hardware. A FASTQ files with 2 billion reads can be processed in about 45 minutes using 8 CPU cores.

# Citation

