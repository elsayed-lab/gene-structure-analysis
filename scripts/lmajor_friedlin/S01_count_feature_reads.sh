#!/usr/bin/zsh

#
# Counts the total number of SL-containing or Poly(A)-containing reads 
#
sl_outfile="../../input/lmajor_hsapiens/S01_num_reads_sl.csv"
polya_outfile="../../input/lmajor_hsapiens/S01_num_reads_polya.csv"
echo "sample_id,num_reads_sl" > $sl_outfile 
echo "sample_id,num_reads_polya" > $polya_outfile 

# SL reads
for x in $SCRATCH/utr_analysis/lmajor-infecting-hsapiens-*/input/HPGL0*; do
	sample_id=$(echo $x | egrep -o "HPGL[0-9]+" | head -1)

    # SL
    infiles=($SCRATCH/utr_analysis/lmajor-infecting-hsapiens-*/*spliced_leader/minlength4_mindiff-2/${sample_id}/results/matched_reads_*.csv)
    num_reads=$(wc -l $infiles | tail -1 | egrep -o "[0-9]+")
	printf "%s,%d\n" $sample_id $num_reads >> $sl_outfile

    # Poly(A)
    infiles=($SCRATCH/utr_analysis/lmajor-infecting-hsapiens-*/poly*/minlength4_mindiff-2/${sample_id}/results/matched_reads_*.csv)
    num_reads=$(wc -l $infiles | tail -1 | egrep -o "[0-9]+")
	printf "%s,%d\n" $sample_id $num_reads >> $polya_outfile
done

