#!/usr/bin/python

import os
import sys
import copy

with open(sys.argv[1], 'r') as in_vcf_file:
    with open('.'.join(sys.argv[1].split('.')[:-1]) + '.lumpy.vcf', 'w') as out_vcf_file:
        line = in_vcf_file.readline().strip()
        header = False
        while line:
            if line.startswith('#'):
                if not header:
                    header_keys = [t.upper() for t in line.split()]
                    seq_key = header_keys.index('SEQ')
                    svlen_key = header_keys.index('SVLEN')
                    svtype_key = header_keys.index('SVTYPE')
                    assert seq_key != -1
                    assert svlen_key != -1
                    assert svtype_key != -1
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
                    out_vcf_file.write('##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description="Secondary breakend in a multi-line variants">\n')
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
                    out_vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')
            else:
                tokens = line.split()
                fields = []
                fields.append(tokens[0])
                fields.append(tokens[1])
                fields.append('.')
                fields.append('.')
                fields.append('<' + tokens[svtype_key] + '>')
                fields.append('30')
                fields.append('PASS')
                info = []
                info.append('END=' + tokens[2])
                info.append('SEQ=' + tokens[seq_key])
                info.append('SVLEN=' + tokens[svlen_key])
                info.append('CIPOS=-10,10')
                info.append('CIEND=-10,10')
                fields.append(';'.join(info) + ';IMPRECISE;SU=0;PE=0;SR=0')
                fields.append('GT:SU:PE:SR')
                #tokens.append('./.:0:0:0')
                out_vcf_file.write('\t'.join(fields) + '\n')
            line = in_vcf_file.readline().strip()
