#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 4 --mem 4gb
PREFIX=Rant
OUTDIR=telomere_reports
mkdir -p $OUTDIR
ls test_spades*/$PREFIX.sorted.fasta | parallel -j 4 python find_telomeres.py {} \> $OUTDIR/{//}.{/.}.telomere_report.txt
