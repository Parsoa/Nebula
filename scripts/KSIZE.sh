export KSZ=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $10 } else if ($0 ~ /.*genotyping.*/ || NF < 8) { print "31k" } else { print $8 } }')
echo KSZ $KSZ
