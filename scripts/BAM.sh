BAM=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.sort.bam" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/reads.bam" } }')
echo BAM $BAM
