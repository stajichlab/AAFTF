#!/bin/bash -l
#SBATCH -N 1 -n 48 -p short --mem 96gb  --out test_spades_nomerge.%A.log
module load AAFTF
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=2
fi
MEM=96
OUTDIR=test_spades_nomerge
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

LEFTTRIMFP=$OUTDIR/${PREFIX}_fastp_1P.fastq.gz
RIGHTTRIMFP=$OUTDIR/${PREFIX}_fastp_2P.fastq.gz

LEFT=$OUTDIR/${PREFIX}_filtered_1.fastq.gz
RIGHT=$OUTDIR/${PREFIX}_filtered_2.fastq.gz


if [ ! -f $LEFT ]; then
    if [ ! -f $LEFTTRIMFP ]; then
	../scripts/AAFTF trim --mem $MEM --method fastp --left $OUTDIR/${SRA}_1.fastq.gz  --right $OUTDIR/${SRA}_2.fastq.gz \
			 -o $OUTDIR/${PREFIX}_fastp -c $CPU --dedup
    fi
    
    ../scripts/AAFTF filter --mem $MEM -c $CPU --left $LEFTTRIMFP --right $RIGHTTRIMFP --aligner bbduk -o $OUTDIR/${PREFIX}
    
    if [ -f $LEFT ]; then
	rm -f $LEFTTRIMFP $RIGHTTRIMFP
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
    ../scripts/AAFTF assemble --mem $MEM --left $LEFT --right $RIGHT -o $ASMFILE -c $CPU -w $OUTDIR/spades --debug
    if [ ! -f $ASMFILE ]; then
	echo "spades failed"
	exit
    fi
fi
if [ ! -f $VECCLEAN ]; then
    ../scripts/AAFTF vecscreen -i $ASMFILE -o $VECCLEAN -c $CPU
fi

if [ ! -f $PURGE ]; then
    ../scripts/AAFTF sourpurge -i $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM --left $LEFT  --right $RIGHT
fi
if [ ! -f $PURGE ]; then
    echo "did not succeed making $PURGE"
    exit
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
