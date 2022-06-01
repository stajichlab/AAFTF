#!/bin/bash -l
#SBATCH -N 1 -n 48 -p short --mem 96gb  --out test_spades_3trim.%A.log
module load AAFTF
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=2
fi
MEM=96
OUTDIR=test_spades_3trim
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

LEFTTRIM=$OUTDIR/${PREFIX}_1P.fastq.gz
RIGHTTRIM=$OUTDIR/${PREFIX}_2P.fastq.gz

LEFTTRIMFP=$OUTDIR/${PREFIX}_fastp2_1P.fastq.gz
RIGHTTRIMFP=$OUTDIR/${PREFIX}_fastp2_2P.fastq.gz

MERGEDTRIM=$OUTDIR/${PREFIX}_fastp_MG.fastq.gz

LEFT=$OUTDIR/${PREFIX}_filtered_1.fastq.gz
RIGHT=$OUTDIR/${PREFIX}_filtered_2.fastq.gz
MERGED=$OUTDIR/${PREFIX}_filtered_U.fastq.gz

if [[ ! -f $LEFT && ! -f $MERGED ]]; then
    if [ ! -f $LEFTTRIMFP ]; then
		../scripts/AAFTF trim --mem $MEM --method fastp \
						 --left $OUTDIR/SRR5223785_1.fastq.gz \
						 --right $OUTDIR/SRR5223785_2.fastq.gz \
						 -o $OUTDIR/${PREFIX}_fastp -c $CPU --dedup --merge

		../scripts/AAFTF trim --mem $MEM --method fastp \
						 --left $OUTDIR/${PREFIX}_fastp_1P.fastq.gz \
						 --right $OUTDIR/${PREFIX}_fastp_2P.fastq.gz \
						 --cutright -o $OUTDIR/${PREFIX}_fastp2 -c $CPU

		../scripts/AAFTF trim --mem $MEM --method bbduk \
						 --left $LEFTTRIMFP \
						 --right $RIGHTTRIMFP \
						 -o $OUTDIR/${PREFIX} -c $CPU
    fi

    ../scripts/AAFTF filter --mem $MEM -c $CPU \
					 --left $LEFTTRIM \
					 --right $RIGHTTRIM \
					 --aligner bbduk -o $OUTDIR/${PREFIX}

    ../scripts/AAFTF filter --mem $MEM -c $CPU \
					 --left $MERGEDTRIM \
					 --aligner bbduk -o $OUTDIR/${PREFIX}

    if [ -f $LEFT ]; then
		rm -f $LEFTTRIM $RIGHTTRIM $MERGEDTRIM $OUTDIR/${PREFIX}_fastp?_[12]P.fastq.gz
    fi
fi

ASMFILE=$OUTDIR/${PREFIX}.spades.fasta
VECCLEAN=$OUTDIR/${PREFIX}.vecscreen.fasta
PURGE=$OUTDIR/${PREFIX}.sourpurge.fasta
CLEANDUP=$OUTDIR/${PREFIX}.rmdup.fasta
PILON=$OUTDIR/${PREFIX}.pilon.fasta
SORTED=$OUTDIR/${PREFIX}.sorted.fasta
STATS=$OUTDIR/${PREFIX}.sorted.stats.txt
if [ ! -f $ASMFILE ]; then
    echo "--mem $MEM --left $LEFT --right $RIGHT --merged $MERGED -o $ASMFILE -c $CPU -w $OUTDIR/spades"
    ../scripts/AAFTF assemble --mem $MEM --left $LEFT --right $RIGHT --merged $MERGED -o $ASMFILE -c $CPU -w $OUTDIR/spades --debug --careful
    if [ ! -f $ASMFILE ]; then
	echo "SPades failed"
	exit
    fi
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

if [ ! -f $PILON ]; then
    ../scripts/AAFTF pilon -i $CLEANDUP -o $PILON -c $CPU --left $LEFT --right $RIGHT --mem $MEM
fi

if [ ! -f $SORTED ]; then
    ../scripts/AAFTF sort -i $CLEANDUP -o $SORTED
fi

if [ ! -f $STATS ] ; then
	../scripts/AAFTF assess -i $SORTED -r $STATS
fi
