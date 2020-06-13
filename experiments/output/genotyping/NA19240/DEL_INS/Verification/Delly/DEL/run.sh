bcftools view genotypes.bcf > genotypes.vcf
../parse_vcf.py genotypes.vcf
tabulate.sh genotypes.bed
