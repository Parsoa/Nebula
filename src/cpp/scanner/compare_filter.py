import sys
import json

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
    with open(f) as json_file:
        return json.load(json_file)

cpp_kmers = read_file(sys.argv[1])
python_kmers = read_file(sys.argv[2])

print(len(cpp_kmers))
print(len(python_kmers))

print('Check CPP vs Python..')
n = 0
for kmer in cpp_kmers:
    if kmer not in python_kmers and reverse_complement(kmer) not in python_kmers:
        print(kmer)
        n += 1
print('Mistakes:', n)

m = 0
print('Check Python vs CPP..')
for kmer in python_kmers:
    if kmer not in cpp_kmers and reverse_complement(kmer) not in cpp_kmers:
        print(kmer)
        m += 1
print('Mistakes:', n)

