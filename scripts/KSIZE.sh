export KSZ=$(echo $P | awk -F/ '{ if ($0 ~ /.*mulation.*/) { print $10 } else { print $8 } }')
echo KSZ $KSZ
