import sys

def canonicalize(seq):
    seq = seq.upper()
    reverse_complement = reverse_complement_sequence(seq)
    return seq if seq < reverse_complement else reverse_complement

def reverse_complement(seq):
    return complement_sequence(seq[::-1])

def reverse_complement_sequence(seq):
    return complement_sequence(seq[::-1])

def complement_sequence(seq):
    # A-> C and C->A
    seq = seq.replace('A', 'Z')
    seq = seq.replace('T', 'A')
    seq = seq.replace('Z', 'T')
    #
    seq = seq.replace('G', 'Z')
    seq = seq.replace('C', 'G')
    seq = seq.replace('Z', 'C')
    #
    return seq

def read_file(f):
    kmers = {}
    with open(f) as counts_file:
        line = counts_file.readline()
        while line:
            tokens = line.split(':')
            kmers[canonicalize(tokens[0])] = int(tokens[1])
            line = counts_file.readline()
    return kmers

cpp_kmers = read_file(sys.argv[1])
python_kmers = read_file(sys.argv[2])

print(len(cpp_kmers))
print(len(python_kmers))

print('Check CPP vs Python..')
n = 0
for kmer in cpp_kmers:
    if cpp_kmers[kmer] != python_kmers[kmer]:
        print(kmer, cpp_kmers[kmer], python_kmers[kmer])
        n += 1

print('Mistakes:', n)
print('Check Python vs CPP..')
m = 0
for kmer in python_kmers:
    if cpp_kmers[kmer] != python_kmers[kmer]:
        m += 1
print('Mistakes:', n)

