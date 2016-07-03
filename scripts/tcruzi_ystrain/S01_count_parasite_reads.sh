#!/usr/bin/env sh
#
# Uses host-filtered fastq reads to count the total number of reads which
# mapped to parasite.
#
infile="../../input/tcruzi_ystrain/S01_num_reads_parasite.csv"
echo "sample_id,num_reads_parasite" > $infile

for x in $SCRATCH/utr_analysis/tcruzi_ystrain/tcruzi-infecting-hsapiens-*/common/HPGL*/tophat/mapped_to_target_untrimmed/align_summary.txt; do
	sample_id=$(echo $x | egrep -o "HPGL[0-9]+" | head -1)
    num_reads=$(tail -4 $x | head -1 | egrep -o "[0-9]+")
	printf "%s,%d\n" $sample_id $num_reads >> $infile
done

