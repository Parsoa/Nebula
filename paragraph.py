import pysam

v = pysam.VariantFile('/share/hormozdiarilab/Codes/NebulousSerendipity/data/HGSV/HG00514.merged_nonredundant.lumpy.short.vcf.gz')
for vcf in v.fetch():
    print(vcf)

