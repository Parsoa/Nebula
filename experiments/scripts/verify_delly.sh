grep -E 'CHROM|DEL' merge.bed > merge.DEL.bed

intersect merge.DEL.bed consensus.DEL.bed > merge.DEL.shared.bed
intersect ./Delly/DEL/genotypes.bed consensus.DEL.bed > delly.DEL.shared.bed

#intersect merge.DEL.bed ./Delly/DEL/genotypes.bed > merge.DEL.shared.bed
#intersect ./Delly/DEL/genotypes.bed merge.DEL.bed -v > delly.DEL.not.bed
#intersect ./Delly/DEL/genotypes.bed merge.DEL.bed > delly.DEL.shared.bed

echo "######################## Deletions ############################"

echo "========================== Delly =============================="
tabulate.sh delly.DEL.shared.bed
verify.sh $1
wc -l *_as_*.bed

echo "========================== Nebula ============================="
tabulate.sh merge.DEL.shared.bed 
verify.sh $1
wc -l *_as_*.bed

echo "======================== Only Delly ==========================="
tabulate.sh delly.DEL.not.bed
verify.sh $1
wc -l *_as_*.bed
