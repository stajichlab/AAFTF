#!/bin/bash

module load AAFTF

OUTDIR=test
PREFIX=Rant
PHYLUM=Ascomycota

mkdir -p $OUTDIR
pushd $OUTDIR
if [ ! -f SRR5223785_1.fastq.gz ]; then
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR522/005/SRR5223785/SRR5223785_1.fastq.gz
fi
if [ ! -f SRR5223785_2.fastq.gz ]; then
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR522/005/SRR5223785/SRR5223785_2.fastq.gz
fi
popd

CPU=24
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
LEFTTRIM=$OUTDIR/${PREFIX}_1P.fastq.gz
RIGHTTRIM=$OUTDIR/${PREFIX}_2P.fastq.gz
LEFT=$OUTDIR/${PREFIX}_filtered_1.fastq.gz
RIGHT=$OUTDIR/${PREFIX}_filtered_2.fastq.gz

if [ ! -f $LEFTTRIM ]; then
	../scripts/AAFTF trim --method bbduk --left $OUTDIR/SRR5223785_1.fastq.gz --right $OUTDIR/SRR5223785_2.fastq.gz -o $OUTDIR/${PREFIX} -c $CPU 
fi
if [ ! -f $LEFT ]; then
 	../scripts/AAFTF filter -c $CPU --left $LEFTTRIM --right $RIGHTTRIM --aligner bbduk -o $OUTDIR/${PREFIX} 
fi
unlink $LEFTTRIM
unlink $RIGHTTRIM

ASMFILE=$OUTDIR/${PREFIX}.spades.fasta
VECCLEAN=$OUTDIR/${PREFIX}.vecscreen.fasta
PURGE=$OUTDIR/${PREFIX}.sourpurge.fasta
CLEANDUP=$OUTDIR/${PREFIX}.rmdup.fasta
PILON=$OUTDIR/${PREFIX}.pilon.fasta
SORTED=$OUTDIR/${PREFIX}.sorted.fasta
STATS=$OUTDIR/${PREFIX}.sorted.stats.txt
if [ ! -f $ASMFILE ]; then
	../scripts/AAFTF assemble --left $LEFT --right $RIGHT -o $ASMFILE -c $CPU --method dipspades
fi
if [ ! -f $VECCLEAN ]; then
	../scripts/AAFTF vecscreen -i $ASMFILE -o $VECCLEAN -c $CPU
fi
if [ ! -f $PURGE ]; then
	../scripts/AAFTF sourpurge -i  $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM --left $LEFT  --right $RIGHT
fi
if [ ! -f $CLEANDUP ]; then
	../scripts/AAFTF rmdup -i $PURGE -o $CLEANDUP -c $CPU -m 1000
fi
if [ ! -f $SORTED ]; then
	../scripts/AAFTF sort -i $CLEANDUP -o $SORTED
fi

if [ ! -f $STATS ] ; then
	../scripts/AAFTF assess -i $SORTED -r $STATS
fi
