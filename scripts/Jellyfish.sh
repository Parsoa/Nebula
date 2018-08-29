JLY=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/mer_counts.jf" } }')
echo JLY $JLY
RJF=$(echo $SIM | awk -v P="$P" -v REF="$REF" -v KSZ="$KSZ" '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/reference_" KSZ ".jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" REF "/mer_counts_" KSZ ".jf" } }')
echo RJF $RJF
