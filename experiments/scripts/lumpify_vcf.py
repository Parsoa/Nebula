#!/usr/bin/python

import os
import sys
import copy

with open(sys.argv[1], 'r') as in_vcf_file:
    with open('.'.join(sys.argv[1].split('.')[:-1]) + '.lumpy.vcf', 'w') as out_vcf_file:
        line = in_vcf_file.readline().strip()
        header = False
        while line:
            if line.startswith('##'):
                if not header:
                    header = True
                    out_vcf_file.write('##fileformat=VCFv4.2\n')
                    out_vcf_file.write('##source=LUMPY\n')
                    out_vcf_file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                    out_vcf_file.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
                    out_vcf_file.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
                    out_vcf_file.write('##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">\n')
                    out_vcf_file.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
                    out_vcf_file.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\n')
                    out_vcf_file.write('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">\n')
                    out_vcf_file.write('##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">\n')
                    out_vcf_file.write('##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">\n')
                    out_vcf_file.write('##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">\n')
                    out_vcf_file.write('##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">\n')
                    out_vcf_file.write('##INFO=<ID=SEQ,Number=1,Type=String,Description="Sequence of the structural variation">\n')
                    out_vcf_file.write('##INFO=<ID=SVCLASS,Number=1,Type=String,Description="Sub-type of the SV">\n')
                    out_vcf_file.write('##INFO=<ID=IS_TRF,Number=1,Type=String,Description="Flag showing if event is in a tandem repeat region.">\n')
                    out_vcf_file.write('##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">\n')
                    out_vcf_file.write('##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">\n')
                    out_vcf_file.write('##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">\n')
                    out_vcf_file.write('##INFO=<ID=BD,Number=.,Type=Integer,Description="Amount of BED evidence supporting the variant across all samples">\n')
                    out_vcf_file.write('##INFO=<ID=EV,Number=.,Type=String,Description="Type of LUMPY evidence contributing to the variant call">\n')
                    out_vcf_file.write('##INFO=<ID=PRPOS,Number=.,Type=String,Description="LUMPY probability curve of the POS breakend">\n')
                    out_vcf_file.write('##INFO=<ID=PREND,Number=.,Type=String,Description="LUMPY probability curve of the END breakend">\n')
                    out_vcf_file.write('##ALT=<ID=DEL,Description="Deletion">\n')
                    out_vcf_file.write('##ALT=<ID=DUP,Description="Duplication">\n')
                    out_vcf_file.write('##ALT=<ID=INV,Description="Inversion">\n')
                    out_vcf_file.write('##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">\n')
                    out_vcf_file.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n')
                    out_vcf_file.write('##ALT=<ID=CNV,Description="Copy number variable region">\n')
                    out_vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    out_vcf_file.write('##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">\n')
                    out_vcf_file.write('##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">\n')
                    out_vcf_file.write('##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">\n')
                    out_vcf_file.write('##FORMAT=<ID=BD,Number=1,Type=Integer,Description="Amount of BED evidence supporting the variant">\n')
            elif line.startswith('#'):
                out_vcf_file.write(line + '\n')
            else:
                tokens = line.split()
                info = tokens[7]
                fields = info.split(';')
                fields.append('CIPOS=-10,10')
                fields.append('CIEND=-10,10')
                keep = []
                has_seq = False
                for i, field in enumerate(fields):
                    key = field.split('=')[0]
                    value = field.split('=')[1]
                    if key == 'SEQ':
                        has_seq = True
                        field = key + '=' + value.upper()
                        if value == '0':
                            line = in_vcf_file.readline().strip()
                            continue
                    if key == 'END':
                        if tokens[4] == '<INS>':
                            field = key + '=' + tokens[1]
                    if key == 'IS_TRF':
                        if value == 'TR':
                            field = key + '=' + 'True'
                        else:
                            field = key + '=' + 'False'
                    if key == 'SVTYPE':
                        svtype = value
                    # discard useless fields
                    if key not in ['CONTIG', 'CONTIG_START', 'CONTIG_END', 'CALLSET', 'UNION']:
                        keep.append(field)
                if svtype != 'NONE':
                    tokens[7] = ';'.join(keep) + ';IMPRECISE;SU=0;PE=0;SR=0'
                    tokens[8] = 'GT:SU:PE:SR'
                    tokens[9] = './.:0:0:0'
                    out_vcf_file.write('\t'.join(tokens) + '\n')
            line = in_vcf_file.readline().strip()
