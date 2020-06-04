grep -E 'CHROM|DEL' merge.bed > merge.DEL.bed
grep -E 'CHROM|INS' merge.bed > merge.INS.bed

intersect merge.DEL.bed consensus.DEL.bed > merge.DEL.shared.bed
intersect ./SVtyper/DEL/genotypes.bed consensus.DEL.bed > svtyper.DEL.shared.bed
intersect svtyper.DEL.shared.bed merge.DEL.bed -v > svtyper.DEL.not.bed

#intersect merge.DEL.bed ./SVtyper/DEL/genotypes.bed > merge.DEL.shared.bed
#intersect ./SVtyper/DEL/genotypes.bed merge.DEL.bed -v > svtyper.DEL.not.bed
#intersect merge.DEL.bed ./SVtyper/DEL/genotypes.bed -b > svtyper.DEL.shared.bed

echo "######################## Deletions ############################"

echo "========================= SVtyper ============================="
tabulate.sh svtyper.DEL.shared.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "========================== Nebula ============================="
tabulate.sh merge.DEL.shared.bed 
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "======================== Only SVtyper ========================="
tabulate.sh svtyper.DEL.not.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed
