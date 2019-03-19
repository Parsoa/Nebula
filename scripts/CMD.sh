CMD=$(echo $PWD | awk '{ if ($0 ~ /.*Simulation.*/) { print "simulate" } else if ($0 ~ /*.genotyping.*/) { print "genotype" } else { print "preprocess" } }')
echo CMD $CMD
