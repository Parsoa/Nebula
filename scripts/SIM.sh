export SIM=$(echo $P   | awk -F/ '{ if ($0 ~ /.*simulation.*/) { print "--simulation " $(NF - 1) } else { print "" } }')
export DESCP=$(echo $P | awk -F/ '{ if ($0 ~ /.*simulation.*/) { print "--description " $(NF - 2) } else { print "" } }')
export DEPTH=$(echo $P | awk -F/ '{ if ($0 ~ /.*simulation.*/) { print "--depth " $(NF - 1) } else { print "" } }')
echo SIM $SIM
echo DESCP $DESCP
echo DEPTH $DEPTH
