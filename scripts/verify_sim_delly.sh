grep -E 'CHROM|DEL' merge.bed > merge.DEL.bed
grep -E 'CHROM|INS' merge.bed > merge.INS.bed

intersect merge.DEL.bed ./Delly/DEL/delly.bed > merge.DEL.shared.bed
#intersect merge.INS.bed ./Delly/INS/delly.bed > merge.INS.shared.bed
intersect ./Delly/DEL/delly.bed merge.DEL.bed -v > delly.DEL.not.bed
#intersect ./Delly/INS/delly.bed merge.INS.bed -v > delly.INS.not.bed
intersect merge.DEL.bed ./Delly/DEL/delly.bed -b > delly.DEL.shared.bed
#intersect merge.INS.bed ./Delly/INS/delly.bed -b > delly.INS.shared.bed

echo "######################## Deletions ############################"

echo "========================== Delly =============================="
tabulate.sh delly.DEL.shared.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "========================== Nebula ============================="
tabulate.sh merge.DEL.shared.bed 
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "======================== Only Delly ==========================="
tabulate.sh delly.DEL.not.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed

#echo "######################## Insertions ############################"
#
#echo "========================= Delly ============================"
#tabulate.sh delly.INS.shared.bed
#verify_sim.sh ../Simulation
#wc -l *_as_*.bed
#
#echo "========================== Nebula =============================="
#tabulate.sh merge.INS.shared.bed 
#verify_sim.sh ../Simulation
#wc -l *_as_*.bed
#
#echo "====================== Only Delly =========================="
#tabulate.sh delly.INS.not.bed
#verify_sim.sh ../Simulation
#wc -l *_as_*.bed
