#! /usr/bin/python

# Converts a VCF file to BAM

import os
import sys
import copy

info_fields = ['END', 'SVLEN', 'SVTYPE']
sample_names = ['HG00514', 'NA19240', 'HG00733']
svtypes = ['INV']

stdout = False
if '-' in sys.argv:
    stdout = True

All = False
if '--all' in sys.argv:
    All = True

repeats = False
if '--repeats' in sys.argv:
    repeats = True

with open(sys.argv[1], 'r') as vcf_file:
    n = 0
    state = 0
    files = {}
    samples = {}
    line = vcf_file.readline()
    while line:
        if state == 0:
            if line.startswith('#CHROM'):
                state = 1
                tokens = line.split()
                for index, name in enumerate(tokens[9:]):
                    if name in sample_names:
                        if not stdout:
                            files[name] = open(os.path.join(os.getcwd(), name + '.INV.bed'), 'w')
                            files[name].write('#CHROM\tBEGIN\t' + '\t'.join([f.upper() for f in info_fields]) + '\tGENOTYPE\n')
                        samples[name] = 9 + index
                if stdout:
                    print('#CHROM\tBEGIN\t' + '\t'.join([f.upper() for f in info_fields]) + '\tGENOTYPE')
        else:
            n += 1
            info = {field.upper(): None for field in info_fields}
            tokens = line.split()
            fields = tokens[7].split(';')
            for field in fields:
                t = field.split('=')
                if t[0] in info:
                    info[t[0]] = t[1]
            if info['SVTYPE'] not in svtypes:
                line = vcf_file.readline()
                #print('SVLEN too short.')
                continue
            if abs(int(info['SVLEN'])) < 10:
                line = vcf_file.readline()
                #print('SVLEN too short.')
                continue
            if all([info[field] for field in info_fields]):
                chrom = 'chr' + tokens[0] if 'chr' not in tokens[0] else tokens[0]
                for sample in samples:
                    vcf_line = chrom + '\t' + tokens[1] + '\t' + '\t'.join([str(info[field]) for field in info_fields])
                    # Fix chrX genotypes, they may look like ./1 sometimes.
                    if tokens[samples[sample]] == '0/0':
                        continue
                    if tokens[samples[sample]] == './.':
                        continue
                    if '1' in tokens[samples[sample]]:
                        if 'x' in chrom.lower():
                            if len(tokens[samples[sample]]) < 3:
                                tokens[samples[sample]] = '1/0'
                    vcf_line += '\t' + tokens[samples[sample]].replace('|', '/') 
                    vcf_line += '\n'
                    if stdout:
                        print(vcf_line.strip())
                    else:
                        files[sample].write(vcf_line)
        line = vcf_file.readline()
