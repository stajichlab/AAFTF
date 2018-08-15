#!/bin/bash

module load ncbi-blast/2.7.1+
module load samtools/1.8
module load bwa/0.7.17
module load trimmomatic
module load java
module load bowtie2
module load python/3
module load SPAdes
module load pilon
module load blobtools

OUTDIR=test
PREFIX=Rant

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

if [ ! -f $OUTDIR/${PREFIX}_1P.fastq ]; then
	./scripts/AAFTF trim --left $OUTDIR/SRR5223785_1.fastq.gz --right $OUTDIR/SRR5223785_2.fastq.gz -o $OUTDIR --prefix $PREFIX \
		-c $CPU --trimmomatic  $TRIMMOMATIC
fi
if [ ! -f $OUTDIR/${PREFIX}_cleaned_1.fq.gz ]; then
 	./scripts/AAFTF filter -i $OUTDIR --prefix $PREFIX --bowtie2 --paired -c $CPU 	
fi
if [ ! -f $OUTDIR/${PREFIX}/scaffolds.fasta ]; then
	./scripts/AAFTF assemble -i $OUTDIR --prefix $PREFIX --paired -c $CPU
fi
if [ ! -f $OUTDIR/${PREFIX}.vecscreen.fasta ]; then
	./scripts/AAFTF vecscreen -i $OUTDIR/$PREFIX/scaffolds.fasta -o $OUTDIR/$PREFIX.vecscreen.fasta
fi
if [ ! -f $OUTDIR/${PREFIX}.vecscreen_purge.fasta ]; then
	./scripts/AAFTF blobpurge -i $OUTDIR/$PREFIX.vecscreen.fasta \
	-o $OUTDIR/$PREFIX.vescreen_purge.fasta --left $OUTDIR/SRR5223785_1.fastq.gz --right $OUTDIR/SRR5223785_2.fastq.gz \
	--phylum Ascomycota 
fi
if [ ! -f $OUTDIR/${PREFIX}.cleaned.fasta ]; then
	./scripts/AAFTF rmdup -i $OUTDIR/$PREFIX.vescreen_purge.fasta -o $OUTDIR/$PREFIX.cleaned.fasta
fi
if [ ! -f $OUTDIR/$PREFIX.cleaned_sorted.fasta ]; then
	./scripts/AAFTF sort -i $OUTDIR/$PREFIX.cleaned.fasta -o $OUTDIR/$PREFIX.cleaned_sorted.fasta
fi
./scripts/AAFTF assess -i $OUTDIR/$PREFIX.vecscreen.fasta -r $OUTDIR/$PREFIX.vecscreen.stats.txt  --stats
#./scripts/AAFTF assess -i $OUTDIR/$PREFIX.vecscreen.fasta -r $OUTDIR/$PREFIX.vecscreen.stats_busco.txt  

