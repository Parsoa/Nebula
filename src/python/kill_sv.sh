# ps -A | grep "break_point"
#for pid in $(ps aux | grep "genotyping" | cut -d " " -f 2); do
#    echo $pid
#    kill -9 $pid
#done
pkill -f "break_point"
pkill -f "prune"
# ps -A | grep "break_point"
