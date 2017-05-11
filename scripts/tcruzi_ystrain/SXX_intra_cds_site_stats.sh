#!/usr/bin/env sh
#
# Processes output files from utr_analysis.py and summarizes statistics for
# SL and Poly(A) sites detected inside annotated CDS regions.
#
wd=$(pwd)

cd /cbcb/nelsayed-scratch/keith/utr_analysis/tcruzi_ystrain

# Created a combined CSV file containing all intra-CDS SL reads
OUTFILE="${wd}/../../input/tcruzi_ystrain/XX_intracds_orfs.csv"

touch $OUTFILE
echo "gene_id,chromosome,strand,position" > $OUTFILE

for stage in tcruzi-*; do
    echo $stage;
    sl_total=0

    for x in ${stage}/spliced_leader/minlength3_mindiff-2/*/results/inside_cds_*; do
        # keep track of number of sites`
        count=$(( $(cat $x | wc -l) - 1 ))
        sl_total=$(( $sl_total + $count ))

        # append to combined output file
        tail -n+2 $x >> $OUTFILE
    done

    echo $sl_total
done

# compress combined output file
gzip $OUTFILE

