RJF=$(echo $SIM | awk -v P="$P" -v REF="$REF" -v KSZ="$KSZ" '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/reference_" KSZ ".jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" REF "/mer_counts.jf" } }')
echo $RJF
