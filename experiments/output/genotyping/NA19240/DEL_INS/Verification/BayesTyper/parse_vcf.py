#!/usr/bin/python

# Converts a VCF file to BAM

import os
import sys
import copy

info_fields = ['END', 'SVTYPE']
svtypes = ['DEL', 'INS', 'INV', 'ALU', 'MEI']

stdout = False
if len(sys.argv) > 2 and sys.argv[2] == '-':
    stdout = True

with open(sys.argv[1], 'r') as vcf_file:
    n = 0
    state = 0
    samples = {}
    line = vcf_file.readline()
    while line:
        if state == 0:
            if line.startswith('#CHROM'):
                state = 1
                tokens = line.split()
                if not stdout:
                    for index, name in enumerate(tokens[9:]):
                        bed_file = open(os.path.join(os.getcwd(), 'genotypes.bed'), 'w')
                        bed_file.write('#CHROM\tBEGIN\t' + '\t'.join([f.upper() for f in info_fields]) + '\tGENOTYPE\n')
                else:
                    print('#CHROM\tBEGIN\t' + '\t'.join([f.upper() for f in info_fields]))
        else:
            n += 1
            tokens = line.split()
            fields = tokens[7].split(';')
            for field in fields:
                t = field.split('=')
                if t[0] == 'VCGR':
                    begin, end = t[1].split('-')
                    begin = begin.split(':')[1]
            svtype = tokens[2][:3]
            chrom = 'chr' + tokens[0] if 'chr' not in tokens[0] else tokens[0]
            if svtype == 'INS':
                end = str(int(end) + 1)
            if not stdout:
                vcf_line = chrom + '\t' + begin + '\t' + end + '\t' + svtype 
                genotype = tokens[9].split(':')[0]
                if '1' in genotype:
                    if 'x' in chrom.lower():
                        if len(genotype) < 3:
                            genotype = '1/0'
                vcf_line += '\t' + genotype 
                vcf_line += '\n'
                bed_file.write(vcf_line)
                #else:
                #    vcf_line = chrom + '\t' + tokens[1] + '\t' + '\t'.join([str(info[field]) for field in info_fields])
                #    print(vcf_line.strip())
        line = vcf_file.readline()
