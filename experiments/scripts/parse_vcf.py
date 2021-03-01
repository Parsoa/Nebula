#! /usr/bin/python

# Converts a VCF file to BAM

import os
import sys
import copy

info_fields = ['END', 'SVTYPE', 'SVLEN', 'SEQ', 'SVCLASS', 'IS_TRF']
svtypes = ['DEL', 'INS', 'INV', 'ALU', 'MEI']

with open(sys.argv[1], 'r') as vcf_file:
    n = 0
    state = 0
    line = vcf_file.readline()
    while line:
        if state == 0:
            if line.startswith('#CHROM'):
                state = 1
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
            # IS_TRF field has three possible values: ., 0, and TR
            info['IS_TRF'] = 'True' if info['IS_TRF'] == 'TR' else 'False'
            if info['IS_TRF'] == 'True':
                if not All and not repeats:
                    line = vcf_file.readline()
                    continue
            else:
                if repeats:
                    line = vcf_file.readline()
                    continue
            if info['SEQ'] == '0':
                line = vcf_file.readline()
                #print('SEQ invalid.')
                continue
            if info['SVTYPE'] not in svtypes:
                line = vcf_file.readline()
                #print('SVTYPE invalid: ', info['SVTYPE'])
                continue
            if abs(int(info['SVLEN'])) < 10:
                line = vcf_file.readline()
                #print('SVLEN too short.')
                continue
            if info['SVTYPE'] == 'DEL':
                info['END'] = int(tokens[1]) + abs(int(info['SVLEN']))
                info['SEQ'] = '.'
            else:
                info['END'] = int(tokens[1]) + 1
            if all([info[field] for field in info_fields]):
                chrom = 'chr' + tokens[0] if 'chr' not in tokens[0] else tokens[0]
                vcf_line = chrom + '\t' + tokens[1] + '\t' + '\t'.join([str(info[field]) for field in info_fields])
                # Fix chrX genotypes, they may look like ./1 sometimes.
                if '1' in tokens[samples[sample]]:
                    if 'x' in chrom.lower():
                        if len(tokens[samples[sample]]) < 3:
                            tokens[samples[sample]] = '1/0'
                vcf_line += '\t' + tokens[samples[sample]].replace('|', '/') 
                vcf_line += '\n'
                print(vcf_line.strip())
        line = vcf_file.readline()
