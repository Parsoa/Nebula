grep -E 'CHROM|DEL' merge.bed > merge.DEL.bed

intersect merge.DEL.bed ./SVtyper/DEL/genotypes.bed > merge.DEL.shared.bed
intersect ./SVtyper/DEL/genotypes.bed merge.DEL.bed -v > svtyper.DEL.not.bed
intersect ./SVtyper/DEL/genotypes.bed merge.DEL.bed > svtyper.DEL.shared.bed

echo "######################## Deletions ############################"

echo "========================= SVtyper ============================="
tabulate.sh svtyper.DEL.shared.bed
verify.sh $1
wc -l *_as_*.bed

echo "========================== Nebula ============================="
tabulate.sh merge.DEL.shared.bed 
verify.sh $1
wc -l *_as_*.bed

echo "======================= Only SVTyper =========================="
tabulate.sh svtyper.DEL.not.bed
verify.sh $1
wc -l *_as_*.bed
