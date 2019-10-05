#!/usr/bin/python

import os
import sys

with open(sys.argv[1]) as vcf_file:
    name = '.'.join(sys.argv[1].split('.')[:-1])
    with open(name + '.bed', 'w') as bed_file:
        line = vcf_file.readline()
        state = 0
        while state == 0:
            if not line.startswith('#'):
                state = 1
                break
            tokens = line.split()
            if 'chrom' in tokens[0].lower():
                state = 1
            line = vcf_file.readline()
        bed_file.write('#CHROM\tBEGIN\tEND\tSVTYPE\tGENOTYPE\n')
        insertions = {}
        while line:
            tokens = line.split()
            info = tokens[7]
            fields = info.split(';')
            svtype = fields[0].split('=')[1]
            if svtype == 'BND':
                if tokens[0] not in insertions:
                    insertions[tokens[0]] = []
                insertions[tokens[0]].append(int(tokens[1]))
                line = vcf_file.readline()
                continue
            for field in fields[1:]:
                key = field.split('=')[0]
                if key == 'END':
                    end = field.split('=')[1]
                    genotype = tokens[9].split(':')[0]
                    tokens[0] = tokens[0] if tokens[0].lower().startswith('chr') else 'chr' + tokens[0]
                    bed_file.write(tokens[0] + '\t' + tokens[1] + '\t' + end + '\t' + svtype + '\t' + genotype + '\n')
                    break
            line = vcf_file.readline()
        found = {}
        for chrom in insertions:
            print(chrom)
            insertions[chrom] = sorted(insertions[chrom])
            i = 0
            while True:
                if i >= len(insertions[chrom]) - 1:
                    break
                s = int(insertions[chrom][i])
                t = int(insertions[chrom][i + 1])
                if int(t - s) < 10000:
                    bed_file.write(chrom + '\t' + str(s) + '\t' + str(t) + '\t' + 'INS' + '\t' + '.\.' + '\n')
                    i += 2
                else:
                    i += 1
