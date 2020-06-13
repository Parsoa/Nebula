#! /usr/bin/python

# Converts a VCF file to BAM

import os
import sys
import copy

info_fields = ['END', 'SVTYPE', 'SEQ']
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
            info = {field.upper(): None for field in info_fields}
            tokens = line.split()
            fields = tokens[7].split(';')
            for field in fields:
                t = field.split('=')
                if t[0] in info:
                    info[t[0]] = t[1]
                elif t[0] == 'CONSENSUS':
                    info['SEQ'] = t[1]
            if info['SVTYPE'] not in svtypes:
                line = vcf_file.readline()
                continue
            if info['SVTYPE'] == 'DEL':
                info['SEQ'] = '.'
                info['SVLEN'] = int(info['END']) - int(tokens[1])
            else:
                info['END'] = int(tokens[1]) + 1
            if all([info[field] for field in info_fields]):
                chrom = 'chr' + tokens[0] if 'chr' not in tokens[0] else tokens[0]
                if not stdout:
                    vcf_line = chrom + '\t' + tokens[1] + '\t' + '\t'.join([str(info[field]) for field in info_fields])
                    tokens[9] = tokens[10].split(':')[0]
                    if '1' in tokens[9]:
                        if 'x' in tokens[0].lower():
                            if len(tokens[9]) < 3:
                                tokens[9] = '1/0'
                    vcf_line += '\t' + tokens[9] 
                    vcf_line += '\n'
                    bed_file.write(vcf_line)
                else:
                    vcf_line = chrom + '\t' + tokens[1] + '\t' + '\t'.join([str(info[field]) for field in info_fields])
                    print(vcf_line.strip())
        line = vcf_file.readline()
