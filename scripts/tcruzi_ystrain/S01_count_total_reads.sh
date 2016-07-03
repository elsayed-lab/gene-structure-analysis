#!/usr/bin/env sh
#
# Counts the total number of reads uses as input for gene structure analysis.
#
outfile="../../input/tcruzi_ystrain/S01_num_reads_total.csv"
echo "sample_id,num_reads_total" > $outfile 

for x in $SCRATCH/utr_analysis/tcruzi_ystrain/*/input/HPGL*/processed/*R1*.gz; do
	sample_id=$(echo $x | egrep -o "HPGL[0-9]+" | head -1)
	num_reads=$(( $(zcat $x | wc -l) / 4 ))
	printf "%s,%d\n" $sample_id $num_reads >> $outfile
done
