#!/bin/bash -l
#SBATCH -N 1 -n 24  --mem 96gb --out test_dipspades.%A.log

module load AAFTF
module load spades/3.12.0
which dipspades.py
which spades.py

MEM=96
OUTDIR=test_dipspades
PREFIX=Rant
PHYLUM=Ascomycota

mkdir -p $OUTDIR

SRA=SRR5223785
FOLDER=$(echo -n $SRA | perl -p -e '$_=sprintf("%s/%03d",substr($_,0,6),substr($_,3,1))')
URL=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${FOLDER}/${SRA}/${SRA}
for DIRECTION in 1 2
do
    if [ ! -f $OUTDIR/${SRA}_${DIRECTION}.fastq.gz ]; then
		echo "downloading $OUTDIR/${SRA}_${DIRECTION}.fastq.gz  from ${URL}_${DIRECTION}.fastq.gz"
		curl -o $OUTDIR/${SRA}_${DIRECTION}.fastq.gz ${URL}_${DIRECTION}.fastq.gz
    fi
done


CPU=24
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
LEFTTRIM=$OUTDIR/${PREFIX}_1P.fastq.gz
RIGHTTRIM=$OUTDIR/${PREFIX}_2P.fastq.gz
MERGEDTRIM=$OUTDIR/${PREFIX}_MG.fastq.gz
LEFTTRIM2=$OUTDIR/${PREFIX}_rnd2_1P.fastq.gz
RIGHTTRIM2=$OUTDIR/${PREFIX}_rnd2_2P.fastq.gz

LEFT=$OUTDIR/${PREFIX}_filtered_1.fastq.gz
RIGHT=$OUTDIR/${PREFIX}_filtered_2.fastq.gz
MERGED=$OUTDIR/${PREFIX}_filtered_U.fastq.gz
if [ ! -f $LEFTTRIM ]; then
    ../scripts/AAFTF trim --method fastp \
		     --left $OUTDIR/SRR5223785_1.fastq.gz \
		     --right $OUTDIR/SRR5223785_2.fastq.gz \
		     -o $OUTDIR/${PREFIX} -c $CPU --merge --dedup
    ../scripts/AAFTF trim --mem $MEM --method bbduk \
		     --left $LEFTTRIM  \
		     --right $RIGHTTRIM  \
		     -o $OUTDIR/${PREFIX}_rnd2 -c $CPU

fi
if [ ! -f $LEFT ]; then
    ../scripts/AAFTF filter --mem $MEM -c $CPU \
		     --left $LEFTTRIM2 --right $RIGHTTRIM2 \
		     --aligner bbduk -o $OUTDIR/${PREFIX}
    ../scripts/AAFTF filter --mem $MEM -c $CPU \
		     --left $MERGEDTRIM \
		     --aligner bbduk -o $OUTDIR/${PREFIX}
fi
if [[ -f $LEFT && -f $RIGHT ]]; then
    rm -f $LEFTTRIM $RIGHTTRIM $MERGEDTRIM $LEFTTRIM2 $RIGHTTRIM2
fi

ASMFILE=$OUTDIR/${PREFIX}.dipspades.fasta
VECCLEAN=$OUTDIR/${PREFIX}.vecscreen.fasta
PURGE=$OUTDIR/${PREFIX}.sourpurge.fasta
CLEANDUP=$OUTDIR/${PREFIX}.rmdup.fasta
PILON=$OUTDIR/${PREFIX}.pilon.fasta
SORTED=$OUTDIR/${PREFIX}.sorted.fasta
STATS=$OUTDIR/${PREFIX}.sorted.stats.txt
if [ ! -f $ASMFILE ]; then
    ../scripts/AAFTF assemble --mem $MEM --left $LEFT --right $RIGHT \
		     --merged $MERGED --out $ASMFILE -c $CPU \
		     --method dipspades --workdir $OUTDIR/dipspades
fi
exit
if [ ! -f $VECCLEAN ]; then
	../scripts/AAFTF vecscreen -i $ASMFILE -o $VECCLEAN -c $CPU
fi
if [ ! -f $PURGE ]; then
    ../scripts/AAFTF sourpurge -i  $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM \
		     --left $LEFT  --right $RIGHT
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
