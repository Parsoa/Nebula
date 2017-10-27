ps -A | grep "genotyping"
#for pid in $(ps aux | grep "genotyping" | cut -d " " -f 2); do
#    echo $pid
#    kill -9 $pid
#done
pkill -f "genotyping"
ps -A | grep "genotyping"
