cat $1 | grep -E 'CHROM|1/1' > v_homozygous.bed
cat $1 | grep -E 'CHROM|1/0|0/1' > v_heterozygous.bed
cat $1 | grep -E 'CHROM|1/0|0/1|1/1' > v_present.bed
cat $1 | grep -E 'CHROM|0/0' > v_00.bed

p="00"
intersect "00.bed" v_present.bed -v -f 0.9 > 00_as_00.bed
intersect "00.bed" v_heterozygous.bed -f 0.9 > 10_as_00.bed
intersect "00.bed" v_homozygous.bed -f 0.9 > 11_as_00.bed
p="10"
intersect "10.bed" v_present.bed -v -f 0.9 > 00_as_10.bed
intersect "10.bed" v_heterozygous.bed -f 0.9 > 10_as_10.bed
intersect "10.bed" v_homozygous.bed -f 0.9 > 11_as_10.bed
p="11"
intersect "11.bed" v_present.bed -v -f 0.9 > 00_as_11.bed
intersect "11.bed" v_heterozygous.bed -f 0.9 > 10_as_11.bed
intersect "11.bed" v_homozygous.bed -f 0.9 > 11_as_11.bed
