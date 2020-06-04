grep -E 'CHROM|DEL' merge.bed > merge.DEL.bed
grep -E 'CHROM|INS' merge.bed > merge.INS.bed

intersect merge.DEL.bed consensus.DEL.bed > merge.DEL.shared.bed
intersect merge.INS.bed consensus.INS.bed > merge.INS.shared.bed
intersect ./BayesTyper/bayestyper_unit_1/genotypes.bed consensus.DEL.bed > bayestyper.DEL.shared.bed
intersect ./BayesTyper/bayestyper_unit_1/genotypes.bed consensus.INS.bed > bayestyper.INS.shared.bed
intersect bayestyper.DEL.shared.bed merge.DEL.bed -v > bayestyper.DEL.not.bed
intersect bayestyper.INS.shared.bed merge.INS.bed -v > bayestyper.INS.not.bed

#intersect merge.DEL.bed ./BayesTyper/DEL/bayestyper.bed > merge.DEL.shared.bed
#intersect merge.INS.bed ./BayesTyper/INS/bayestyper.bed > merge.INS.shared.bed
#intersect ./BayesTyper/DEL/bayestyper.bed merge.DEL.bed -v > bayestyper.DEL.not.bed
#intersect ./BayesTyper/INS/bayestyper.bed merge.INS.bed -v > bayestyper.INS.not.bed
#intersect merge.DEL.bed ./BayesTyper/DEL/bayestyper.bed -b > bayestyper.DEL.shared.bed
#intersect merge.INS.bed ./BayesTyper/INS/bayestyper.bed -b > bayestyper.INS.shared.bed

echo "######################## Deletions ############################"

echo "======================== BayesTyper ============================"
tabulate.sh bayestyper.DEL.shared.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "========================== Nebula =============================="
tabulate.sh merge.DEL.shared.bed 
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "====================== Only BayesTyper =========================="
tabulate.sh bayestyper.DEL.not.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "######################## Insertions ############################"

echo "========================= BayesTyper ============================"
tabulate.sh bayestyper.INS.shared.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "========================== Nebula =============================="
tabulate.sh merge.INS.shared.bed 
verify_sim.sh ../Simulation
wc -l *_as_*.bed

echo "====================== Only BayesTyper =========================="
tabulate.sh bayestyper.INS.not.bed
verify_sim.sh ../Simulation
wc -l *_as_*.bed
