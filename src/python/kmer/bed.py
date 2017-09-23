import pybedtools
from Bio import SeqIO

def read_tracks_from_bed_file(file) :
    bedtool = pybedtools.BedTool(file)
    for track in bedtool :
        print('chrom:', track.chrom, '[', '{:12d}'.format(int(track.start)), ',', '{:12d}'.format(int(track.end)), ']')
        extract_reference_sequence(track)

def extract_reference_sequence(track) :
    fasta = pybedtools.BedTool.sequence('/Users/parsoakhorsand/Davis/Projects/NebulousSerendipity/data/hg38.fa')
    sequence = track.sequence(fi = fasta)
    print(open(track.seqfn).read())
    # Change this later
    # fasta = SeqIO.parse(open('/Users/parsoakhorsand/Davis/Projects/NebulousSerendipity/data/hg38.fa'), 'fasta')
    # for read in fasta :
    #     name, sequence = read.id, read.seq.tostring()
