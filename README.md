### AAFTF - Automatic Assembly For The Fungi
*Jason Stajich and Jon Palmer*

Requirements
===================
- BBTools - https://jgi.doe.gov/data-and-tools/bbtools/
- Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic (Optional)
- bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml (Optional)
- bwa - https://github.com/lh3/bwa
- Pilon - https://github.com/broadinstitute/pilon/wiki
- sourmash - https://sourmash.readthedocs.io/ (install via conda/pip)
- NCBI BLAST+ - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
- minimap2 - https://github.com/lh3/minimap2
Assemblers
- SPAdes - http://cab.spbu.ru/software/spades/
- megahit - https://github.com/voutcn/megahit
- dipspades - (SPAdes 3.11.1 - note it is not part of later SPAdes packages) http://cab.spbu.ru/files/release3.11.1/dipspades_manual.html

Authors
============
* Jason Stajich [@hyphaltip](https://github.com/hyphaltip) - http://lab.stajich.org
* Jon Palmer [@nextgenusfs](https://github.com/nextgenusfs) - https://twitter.com/jonpalmer2013

Notes
===========
This is partially a python re-write of [JAAWS](https://github.com/nextgenusfs/jaaws) which was a unix shell based cleanup and assembly tool written by Jon.

Steps / Procedures
==================
1. trim                Trim FASTQ input reads - with BBMap
2. filter              Filter contaminanting reads - with BBMap
3. assemble            Assemble reads - with SPAdes
4. vecscreen           Vector and Contaminant Screening of assembled contigs - with BlastN based method to replicate NCBI screening
5. sourpurge           Purge contigs based on sourmash results - with sourmash
6. rmdup               Remove duplicate contigs - using minimap2 to find duplicates
7. pilon               Polish contig sequences with Pilon - uses Pilon
8. sort                Sort contigs by length and rename FASTA headers
9. assess              Assess completeness of genome assembly
10. pipeline           Run AAFTF pipeline all in one go.


# Typical runs


## Trimming and Filtering

```
MEM=128 # 128gb
BASE=STRAINX
READSDIR=reads
TRIMREAD=reads_trimmed
CPU=8
AAFTF trim --method bbduk --memory $MEM -c $CPU \
 --left $READSDIR/${BASE}_R1.fq.gz --right $READSDIR/${BASE}_R2.fq.gz \
  -o $TRIMREAD/${BASE}
# this step make take a lot of memory depending on how many filtering libraries you use
AAFTF filter -c $CPU --memory $MEM --aligner bbduk \
	  -o $TRIMREAD/${BASE} --left $TRIMREAD/${BASE}_1P.fastq.gz --right $TRIMREAD/${BASE}_2P.fastq.gz
```

## Assembly

```
CPU=24
MEM=96
LEFT=$TRIMREAD/${BASE}_filtered_1.fastq.gz
RIGHT=$TRIMREAD/${BASE}_filtered_2.fastq.gz
WORKDIR=working_AAFTF
OUTDIR=genomes
ASMFILE=$OUTDIR/${BASE}.spades.fasta
mkdir -p $WORKDIR $OUTDIR
AAFTF assemble -c $CPU --mem $MEM \
	  --left $LEFT --right $RIGHT  \
	   -o $ASMFILE -w $WORKDIR/spades_$BASE
```

## vectrim

```
CPU=16
MEM=16
LEFT=$TRIMREAD/${BASE}_filtered_1.fastq.gz
RIGHT=$TRIMREAD/${BASE}_filtered_2.fastq.gz
WORKDIR=working_AAFTF
OUTDIR=genomes
ASMFILE=$OUTDIR/${BASE}.spades.fasta
VECTRIM=$OUTDIR/${BASE}.vecscreen.fasta
mkdir -p $WORKDIR $OUTDIR
AAFTF vecscreen -c $CPU -i $ASMFILE -o $VECTRIM
```
