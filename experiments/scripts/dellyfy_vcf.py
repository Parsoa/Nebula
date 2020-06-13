#!/usr/bin/python

import os
import sys
import copy

id_counters = {}

with open(sys.argv[1], 'r') as in_vcf_file:
    with open('.'.join(sys.argv[1].split('.')[:-1]) + '.delly.vcf', 'w') as out_vcf_file:
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
                    #out_vcf_file.write('##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">\n')
                    out_vcf_file.write('##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">\n')
                    #out_vcf_file.write('##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">\n')
                    #out_vcf_file.write('##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">\n')
                    out_vcf_file.write('##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">\n')
                    out_vcf_file.write('##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">\n')
                    out_vcf_file.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
                    out_vcf_file.write('##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">\n')
                    out_vcf_file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                    #out_vcf_file.write('##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">\n')
                    out_vcf_file.write('##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">\n')
                    out_vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    #out_vcf_file.write('##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled genotype likelihoods for RR,RA,AA genotypes">\n')
                    #out_vcf_file.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
                    #out_vcf_file.write('##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">\n')
                    #out_vcf_file.write('##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Raw high-quality read counts for the SV">\n')
                    #out_vcf_file.write('##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Raw high-quality read counts for the left control region">\n')
                    #out_vcf_file.write('##FORMAT=<ID=RCR,Number=1,Type=Integer,Description="Raw high-quality read counts for the right control region">\n')
                    #out_vcf_file.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Read-depth based copy-number estimate for autosomal sites">\n')
                    #out_vcf_file.write('##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">\n')
                    #out_vcf_file.write('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">\n')
                    #out_vcf_file.write('##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">\n')
                    #out_vcf_file.write('##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">\n')
                    for c in range(1, 23):
                        out_vcf_file.write('##contig=<ID=chr' + str(c) + '>\n')
                    out_vcf_file.write('##contig=<ID=chrX>\n')
            elif line.startswith('#'):
                tokens = line.split()[:10]
                out_vcf_file.write('\t'.join(tokens[:10]) + '\n')
            else:
                tokens = line.split()[:10]
                info = tokens[7]
                fields = info.split(';')
                fields.append('CIPOS=-10,10')
                fields.append('CIEND=-10,10')
                keep = []
                has_seq = False
                for i, field in enumerate(fields):
                    if not '=' in field:
                        continue
                    key = field.split('=')[0]
                    value = field.split('=')[1]
                    if key == 'SEQ':
                        has_seq = True
                        field = 'CONSENSUS' + '=' + value.upper()
                        if value == '0':
                            line = in_vcf_file.readline().strip()
                            continue
                    if key == 'END':
                        if tokens[4] == '<INS>':
                            field = key + '=' + tokens[1]
                    if key == 'SVTYPE':
                        svtype = value
                        if svtype not in id_counters:
                            id_counters[svtype] = 1000
                        else:
                            id_counters[svtype] += 1
                    # discard useless fields
                    if key not in ['CONTIG', 'CONTIG_START', 'CONTIG_END', 'CALLSET', 'UNION', 'SVCLASS', 'SVLEN', 'IS_TRF', 'SU', 'PE', 'SR']:
                        keep.append(field)
                if svtype != 'NONE':
                    tokens[2] = svtype + '{:08d}'.format(id_counters[svtype])
                    if svtype == 'DEL':
                        tokens[3] = 'N'
                    tokens[6] = '.'
                    tokens[7] = ';'.join(keep) + ';IMPRECISE;MAPQ=60;CT=3to5;INSLEN=0;CHR2=' + tokens[0]
                    tokens[8] = 'GT'
                    tokens[9] = './.'
                    out_vcf_file.write('\t'.join(tokens) + '\n')
            line = in_vcf_file.readline().strip()
