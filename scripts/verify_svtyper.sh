grep SVTYPE=DEL svtyper.vcf > svtyper.DEL.vcf

parse_svtyper_ouput.py svtyper.DEL.vcf

grep DEL merge.bed > merge.DEL.bed

echo 'Intersecting..'
intersect merge.DEL.bed svtyper.DEL.bed > merge.DEL.shared.bed
intersect merge.DEL.bed svtyper.DEL.bed -b > svtyper.DEL.shared.bed

echo "SVtyper results: =========================================="
tabulate_svtyper.sh svtyper.DEL.shared.bed
verify.sh $1
wc -l *.bed

echo "Nebula results: =========================================="
tabulate.sh merge.DEL.shared.bed 
verify.sh $1
wc -l *.bed
