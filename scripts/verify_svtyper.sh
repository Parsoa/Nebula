grep DEL merge.bed > merge.DEL.bed

intersect merge.DEL.bed svtyper.DEL.bed -i > merge.DEL.shared.bed
intersect merge.DEL.bed svtyper.DEL.bed -i -b > svtyper.shared.bed

echo "SVtyper results: =========================================="
svtyper_tabulate.sh svtyper.shared.bed
verify.sh $1
wc -l *.bed

echo "Nebula results: =========================================="
tabulate.sh merge.DEL.shared.bed 
verify.sh $1
wc -l *.bed

