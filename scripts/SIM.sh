export SIM=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--simulation " $9 } else { print "" } }')
export DESCP=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--description " $8 } else { print "" } }')
export DEPTH=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print "--depth " $9 } else { print "" } }')
echo SIM $SIM
echo DESCP $DESCP
echo DEPTH $DEPTH
