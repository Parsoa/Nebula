import khmer

if __name__ == '__main__':
    counttable = khmer.Counttable(31, 16e9, 4)
    n, kmers = counttable.consume_seqfile('/share/hormozdiarilab/Data/Genomes/Illumina/1KG_Trio/HG00513.fq')
