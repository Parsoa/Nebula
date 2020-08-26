cat $1 | grep -E 'CHROM|1/1'         > v_homozygous.bed
cat $1 | grep -E 'CHROM|1/0|0/1'     > v_heterozygous.bed
cat $1 | grep -E 'CHROM|1/0|0/1|1/1' > v_present.bed
cat $1 | grep -E 'CHROM|0/0'         > v_absent.bed
cat $1 | grep -E 'CHROM|\./\.'       > v_na.bed

p="00"
intersect "00.bed" v_na.bed > NA_as_00.bed
intersect "00.bed" v_present.bed -v > 00_as_00.bed
intersect "00.bed" v_heterozygous.bed > 10_as_00.bed
intersect "00.bed" v_homozygous.bed > 11_as_00.bed
p="10"
intersect "10.bed" v_na.bed > NA_as_10.bed
intersect "10.bed" v_present.bed -v > 00_as_10.bed
intersect "10.bed" v_heterozygous.bed > 10_as_10.bed
intersect "10.bed" v_homozygous.bed > 11_as_10.bed
p="11"
intersect "11.bed" v_na.bed > NA_as_11.bed
intersect "11.bed" v_present.bed -v > 00_as_11.bed
intersect "11.bed" v_heterozygous.bed > 10_as_11.bed
intersect "11.bed" v_homozygous.bed > 11_as_11.bed

#cat $1 | awk '{if ($1 == "#CHROM" || $NF == "0/0") {print $0} }' > v_absent.bed
#cat $1 | awk '{if ($1 == "#CHROM" || $NF == "1/0" || $NF == "0/1") {print $0} }' > v_heterozygous.bed
#cat $1 | awk '{if ($1 == "#CHROM" || $NF == "1/1") {print $0} }' > v_homozygous.bed
#cat $1 | awk '{if ($1 == "#CHROM" || $NF == "./.") {print $0} }' > v_na.bed
#
#p="00"
#intersect "00.bed" v_na.bed > NA_as_00.bed
#intersect "00.bed" v_absent.bed > 00_as_00.bed
#intersect "00.bed" v_heterozygous.bed > 10_as_00.bed
#intersect "00.bed" v_homozygous.bed > 11_as_00.bed
#p="10"
#intersect "10.bed" v_na.bed > NA_as_10.bed
#intersect "10.bed" v_absent.bed > 00_as_10.bed
#intersect "10.bed" v_heterozygous.bed > 10_as_10.bed
#intersect "10.bed" v_homozygous.bed > 11_as_10.bed
#p="11"
#intersect "11.bed" v_na.bed > NA_as_11.bed
#intersect "11.bed" v_absent.bed > 00_as_11.bed
#intersect "11.bed" v_heterozygous.bed > 10_as_11.bed
#intersect "11.bed" v_homozygous.bed > 11_as_11.bed
