(head -n 1 Paragraph/DEL/present.bed && tail -q -n +2 ./Paragraph/DEL/present.bed ./Delly/DEL/present.bed ./SVtyper/DEL/present.bed | sort -u -k1,1V -k2,2n) > consensus.DEL.bed
(head -n 1 Paragraph/INS/present.bed && tail -q -n +2 ./Paragraph/INS/present.bed |sort -u -k1,1V -k2,2n) > consensus.INS.bed
sed -i '/^[[:space:]]*$/d' consensus.DEL.bed
sed -i '/^[[:space:]]*$/d' consensus.INS.bed
