ps -A | grep "python -m kmer.sv"
# for pid in $(ps -A | grep "python -m kmer.count_server" | cut -c 1-5); do
#     kill -9 $pid
# done
pkill -f "count_server"
ps -A | grep "python -m kmer.sv"