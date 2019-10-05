grep DEL merge.bed > merge.DEL.bed

intersect merge.DEL.bed svtyper.DEL.bed > merge.DEL.shared.bed
intersect merge.DEL.bed svtyper.DEL.bed -b > svtyper.shared.bed

echo "SVtyper results: =========================================="
svtyper_tabulate.sh svtyper.shared.bed
verify_sim.sh $1
wc -l *.bed

echo "Nebula results: =========================================="
tabulate.sh merge.DEL.shared.bed 
verify_sim.sh $1
wc -l *.bed

