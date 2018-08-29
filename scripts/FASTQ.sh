FSQ=$(echo $SIM | awk -v P=$P -v GEN=$GEN '{ if ($0 ~ /.*simulation.*/) { print P "/../Simulation/test.fq" } else { print "/share/hormozdiarilab/Experiments/Jellyfish/" GEN "/reads.fq" } }')
echo FSQ $FSQ
