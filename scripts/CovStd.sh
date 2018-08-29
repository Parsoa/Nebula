DEP=$(echo $P | awk -v P=$P '{ if ($0 ~ /.*simulation.*/) { print P "/../UniqueKmersDepthOfCoverageEstimationJob/stats.json"} else { print P "/../../UniqueKmersDepthOfCoverageEstimationJob/stats.json" } }')
echo $DEP
COV="$(head -n 10 "$DEP" | sed -n 's/\s*"mean":\s*\([0-9]*\.[0-9]*\).*/\1/p')"
STD="$(head -n 10 "$DEP" | sed -n 's/\s*"std":\s*\([0-9]*\.[0-9]*\).*/\1/p')"
echo $COV
echo $STD
