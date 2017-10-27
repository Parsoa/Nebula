ps -A | grep "python -m kmer.sv"
for pid in $(ps -A | grep "python -m kmer.genotyping" | cut -c 1-5); do
    kill -9 $pid
done
ps -A | grep "python -m kmer.sv"
