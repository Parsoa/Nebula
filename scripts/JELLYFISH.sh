RJF=$(echo $SIM | awk -v P="$P" -v REF="$REF" '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/reference_32k.jf" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" REF "/mer_counts_32k.jf" } }')
RJF="/share/hormozdiarilab/Experiments/Jellyfish/hg19/mer_counts_32k.jf"
echo RJF $RJF
