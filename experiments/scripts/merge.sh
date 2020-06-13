#!/bin/bash
module load bcftools/1.10.2
(head -n 44 HG00514.merged_nonredundant.unified.paragraph.vcf && tail -n +45 HG00514.merged_nonredundant.unified.paragraph.vcf | sort -k1,1V -k2,2n) > HG00514.merged_nonredundant.unified.paragraph.sorted.vcf
(head -n 44 HG00733.merged_nonredundant.unified.paragraph.vcf && tail -n +45 HG00733.merged_nonredundant.unified.paragraph.vcf | sort -k1,1V -k2,2n) > HG00733.merged_nonredundant.unified.paragraph.sorted.vcf

bgzip -c --index HG00514.merged_nonredundant.unified.paragraph.sorted.vcf > HG00514.merged_nonredundant.unified.paragraph.sorted.vcf.gz
bgzip -c --index HG00733.merged_nonredundant.unified.paragraph.sorted.vcf > HG00733.merged_nonredundant.unified.paragraph.sorted.vcf.gz

bcftools index HG00514.merged_nonredundant.unified.paragraph.sorted.vcf.gz
bcftools index HG00733.merged_nonredundant.unified.paragraph.sorted.vcf.gz
bcftools merge --force-samples HG00514.merged_nonredundant.unified.paragraph.sorted.vcf.gz HG00733.merged_nonredundant.unified.paragraph.sorted.vcf.gz > HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.vcf

grep -E "#|DEL" HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.vcf > HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.DEL.vcf
grep -E "#|INS" HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.vcf > HG00514_HG00733.merged_nonredundant.unified.paragraph.sorted.INS.vcf
