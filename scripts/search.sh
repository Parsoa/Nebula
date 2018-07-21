#!/bin/bash
P=$(pwd)
echo $P
cd $P
for c in $(seq 1 22); do
    echo chr"$c"
    grep --line-buffered $1 chr"$c"_strand_"$i".fa | awk -v c=$c '{print c "     " $0}'
    for i in $(seq 1 2); do
        echo strand"$i"
        grep --line-buffered $1 chr"$c"_strand_"$i".fa | awk -v c=$c '{print c "     " $0}'
    done
done
exit
for c in $(seq 5 8); do
    for i in $(seq 1 2); do
        for j in $(seq 1 2); do
            grep --line-buffered $1 chr"$c"_strand_"$i"."$j".fq | awk -v c=$c '{print c "     " $0}' &
        done
    done
done
for c in $(seq 9 12); do
    for i in $(seq 1 2); do
        for j in $(seq 1 2); do
            grep --line-buffered $1 chr"$c"_strand_"$i"."$j".fq | awk -v c=$c '{print c "     " $0}' &
        done
    done
done
for c in $(seq 13 16); do
    for i in $(seq 1 2); do
        for j in $(seq 1 2); do
            grep --line-buffered $1 chr"$c"_strand_"$i"."$j".fq | awk -v c=$c '{print c "     " $0}' &
        done
    done
done
for c in $(seq 17 20); do
    for i in $(seq 1 2); do
        for j in $(seq 1 2); do
            grep --line-buffered $1 chr"$c"_strand_"$i"."$j".fq | awk -v c=$c '{print c "     " $0}' &
        done
    done
done
for c in $(seq 21 22); do
    for i in $(seq 1 2); do
        for j in $(seq 1 2); do
            grep --line-buffered $1 chr"$c"_strand_"$i"."$j".fq | awk -v c=$c '{print c "     " $0}' &
        done
    done
done
