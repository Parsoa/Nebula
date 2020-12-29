grep -E 'CHROM|DEL' merge.bed > merge.DEL.bed
grep -E 'CHROM|INS' merge.bed > merge.INS.bed

intersect merge.DEL.bed consensus.DEL.bed > merge.DEL.shared.bed
intersect merge.INS.bed consensus.INS.bed > merge.INS.shared.bed
intersect merge.INS.bed consensus.MEI.bed > merge.INS.shared.bed
intersect consensus.DEL.bed ./BayesTyper/bayestyper_unit_1/genotypes.bed -b > bayestyper.DEL.shared.bed
intersect consensus.INS.bed ./BayesTyper/bayestyper_unit_1/genotypes.bed -b > bayestyper.INS.shared.bed
intersect consensus.MEI.bed ./BayesTyper/bayestyper_unit_1/genotypes.bed -b > bayestyper.INS.shared.bed

#intersect merge.DEL.bed ./BayesTyper/bayestyper_unit_1/genotypes.bed > merge.DEL.shared.bed
#intersect merge.INS.bed ./BayesTyper/bayestyper_unit_1/genotypes.bed > merge.INS.shared.bed
#intersect ./BayesTyper/bayestyper_unit_1/genotypes.bed merge.DEL.bed > bayestyper.DEL.shared.bed
#intersect ./BayesTyper/bayestyper_unit_1/genotypes.bed merge.INS.bed > bayestyper.INS.shared.bed
#intersect ./Paragraph/DEL/genotypes.bed merge.DEL.bed -v > bayestyper.DEL.not.bed
#intersect ./Paragraph/INS/genotypes.bed merge.INS.bed -v > bayestyper.INS.not.bed

echo "######################## Deletions #############################"

echo "======================== BayesTyper ============================"
tabulate.sh bayestyper.DEL.shared.bed
verify.sh $1
wc -l *_as_*.bed

echo "========================== Nebula =============================="
tabulate.sh merge.DEL.shared.bed 
verify.sh $1
wc -l *_as_*.bed

echo "######################## Insertions ############################"

echo "========================= BayesTyper ==========================="
tabulate.sh bayestyper.INS.shared.bed
verify.sh $1
wc -l *_as_*.bed

echo "========================== Nebula =============================="
tabulate.sh merge.INS.shared.bed 
verify.sh $1
wc -l *_as_*.bed
