#!/usr/bin/env sh
#
# Counts the total number of reads uses as input for gene structure analysis.
#
outfile="../../input/lmajor_hsapiens/S01_num_reads_total.csv"
echo "sample_id,num_reads_total" > $outfile 

for x in $SCRATCH/utr_analysis/lmajor-infecting-hsapiens-*/input/HPGL*/processed/*R1*.gz; do
	sample_id=$(echo $x | egrep -o "HPGL[0-9]+" | head -1)
	num_reads_r1=$(( $(zcat $x | wc -l) / 4 ))
	num_reads_r2=$(( $(zcat ${x/R1/R2} | wc -l) / 4 ))
    printf "%s,%d\n" $sample_id $(($num_reads_r1 + $num_reads_r2)) >> $outfile
done
