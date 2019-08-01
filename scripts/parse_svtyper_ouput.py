#!/usr/bin/python

import os
import sys

with open(sys.argv[1]) as vcf_file:
    with open(sys.argv[1] + '.bed', 'w') as bed_file:
        line = vcf_file.readline()
        state = 0
        while state == 0:
            tokens = line.split()
            if 'chrom' in tokens[0].lower():
                state = 1
            line = vcf_file.readline()
        bed_file.write('CHROM\tBEGIN\tEND\tSVTYPE\tGENOTYPE\n')
        while line:
            tokens = line.split()
            info = tokens[7]
            fields = info.split(';')
            svtype = fields[0].split('=')[1]
            for field in fields[1:]:
                key = field.split('=')[0]
                if key == 'END':
                    end = field.split('=')[1]
                    genotype = tokens[9].split(':')[0]
                    bed_file.write('chr' + tokens[0] + '\t' + tokens[1] + '\t' + end + '\t' + svtype + '\t' + genotype + '\n')
                    break
            line = vcf_file.readline()
