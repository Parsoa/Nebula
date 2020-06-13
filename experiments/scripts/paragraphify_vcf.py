#!/usr/bin/python

import os
import re
import sys
import copy

id_counters = {}

with open(sys.argv[1], 'r') as in_vcf_file:
    with open('.'.join(sys.argv[1].split('.')[:-1]) + '.paragraph.vcf', 'w') as out_vcf_file:
        line = in_vcf_file.readline().strip()
        header = False
        while line:
            if line.startswith('##'):
                if not header:
                    header = True
                    out_vcf_file.write('##fileformat=VCFv4.1\n')
                    out_vcf_file.write('##fileDate=20200419\n')
                    out_vcf_file.write('##ALT=<ID=DEL,Description="Deletion">\n')
                    out_vcf_file.write('##ALT=<ID=DUP,Description="Duplication">\n')
                    out_vcf_file.write('##ALT=<ID=INV,Description="Inversion">\n')
                    out_vcf_file.write('##ALT=<ID=TRA,Description="Translocation">\n')
                    out_vcf_file.write('##ALT=<ID=INS,Description="Insertion">\n')
                    out_vcf_file.write('##FILTER=<ID=LowQual,Description="PE/SR support below 3 or mapping quality below 20.">\n')
                    out_vcf_file.write('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">\n')
                    out_vcf_file.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">\n')
                    out_vcf_file.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">\n')
                    out_vcf_file.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n')
                    out_vcf_file.write('##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">\n')
                    out_vcf_file.write('##INFO=<ID=SEQ,Number=1,Type=String,Description="Split-read consensus sequence">\n')
                    out_vcf_file.write('##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">\n')
                    out_vcf_file.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
                    out_vcf_file.write('##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">\n')
                    out_vcf_file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                    out_vcf_file.write('##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">\n')
                    out_vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    for c in range(1, 23):
                        out_vcf_file.write('##contig=<ID=chr' + str(c) + '>\n')
                    out_vcf_file.write('##contig=<ID=chrX>\n')
            elif line.startswith('#'):
                tokens = line.split()[:10]
                out_vcf_file.write('\t'.join(tokens[:10]) + '\n')
            else:
                tokens = line.split()[:10]
                info = tokens[7]
                fields = []
                has_seq = False
                for i, field in enumerate(info.split(';')):
                    if not '=' in field:
                        continue
                    key = field.split('=')[0]
                    value = field.split('=')[1]
                    if key == 'SEQ':
                        if re.match('[ATCGN]+', value.upper()):
                            has_seq = True
                            fields.append(field)
                    if key == 'END':
                        if tokens[4] == '<INS>':
                            field = key + '=' + tokens[1]
                        fields.append(field)
                    if key == 'SVTYPE':
                        svtype = value
                        if svtype not in id_counters:
                            id_counters[svtype] = 1000
                        else:
                            id_counters[svtype] += 1
                        fields.append(field)
                fields.append('CIPOS=-10,10')
                fields.append('CIEND=-10,10')
                if not has_seq:
                    if svtype == 'INS':
                        line = in_vcf_file.readline().strip()
                        continue
                if svtype != 'NONE':
                    tokens[2] = svtype + '{:08d}'.format(id_counters[svtype])
                    tokens[5] = '30'
                    tokens[6] = '.'
                    tokens[7] = ';'.join(fields) + ';IMPRECISE;MAPQ=60'
                    tokens[8] = 'GT'
                    tokens[9] = './.'
                    out_vcf_file.write('\t'.join(tokens) + '\n')
            line = in_vcf_file.readline().strip()
