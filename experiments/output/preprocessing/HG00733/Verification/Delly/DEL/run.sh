bcftools view genotypes.bcf > genotypes.vcf
../parse_vcf.py genotypes.vcf
tabulate.sh genotypes.bed
verify.sh HG00514.unified.bed
