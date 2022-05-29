#!/bin/bash -l
#SBATCH -N 1 -n 2 -p short --mem 24gb  --out test_mito.%A.log
module load AAFTF
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=2
fi
MEM=24
OUTDIR=test_mito
PREFIX=Rant
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

MITO=$OUTDIR/${PREFIX}.mito.fasta
STATS=$OUTDIR/${PREFIX}.mito.stats.txt
if [ ! -f $MITO ]; then
    ../scripts/AAFTF mito --left $LEFT --right $RIGHT -o $MITO
fi

if [[ -f $MITO && ! -f $STATS ]]; then
    ../scripts/AAFTF assess -i $MITO -r $STATS
fi
