### AAFTF - Automatic Assembly For The Fungi
*Jason Stajich and Jon Palmer*

![AAFTF logo](docs/AAFTF.png)

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

Trimming options spelled out:
```
usage: AAFTF trim [-h] [-q] [-o BASENAME] [-c cpus] [-ml MINLEN] -l LEFT
                  [-r RIGHT] [-v] [--pipe] [--method {bbduk,trimmomatic}]
                  [-m MEMORY] [--trimmomatic trimmomatic_jar]
                  [--trimmomatic_adaptors TRIMMOMATIC_ADAPTORS]
                  [--trimmomatic_clip TRIMMOMATIC_CLIP]
                  [--trimmomatic_leadingwindow TRIMMOMATIC_LEADINGWINDOW]
                  [--trimmomatic_trailingwindow TRIMMOMATIC_TRAILINGWINDOW]
                  [--trimmomatic_slidingwindow TRIMMOMATIC_SLIDINGWINDOW]
                  [--trimmomatic_quality TRIMMOMATIC_QUALITY]

This comamnd trims reads in FASTQ format to remove low quality reads and trim
adaptor sequences

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           Do not output warnings to stderr
  -o BASENAME, --out BASENAME
                        Output basename, default to base name of --left reads
  -c cpus, --cpus cpus  Number of CPUs/threads to use.
  -ml MINLEN, --minlen MINLEN
                        Minimum read length after trimming, default: 75
  -l LEFT, --left LEFT  left/forward reads of paired-end FASTQ or single-end
                        FASTQ.
  -r RIGHT, --right RIGHT
                        right/reverse reads of paired-end FASTQ.
  -v, --debug           Provide debugging messages
  --pipe                AAFTF is running in pipeline mode
  --method {bbduk,trimmomatic}
                        Program to use for adapter trimming
  -m MEMORY, --memory MEMORY
                        Max Memory (in GB)
  --trimmomatic trimmomatic_jar, --jar trimmomatic_jar
                        Trimmomatic JAR path

Trimmomatic options:
  Trimmomatic trimming options

  --trimmomatic_adaptors TRIMMOMATIC_ADAPTORS
                        Trimmomatic adaptor file, default: TruSeq3-PE.fa
  --trimmomatic_clip TRIMMOMATIC_CLIP
                        Trimmomatic clipping, default:
                        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
  --trimmomatic_leadingwindow TRIMMOMATIC_LEADINGWINDOW
                        Trimmomatic window processing arguments, default:
                        LEADING:3
  --trimmomatic_trailingwindow TRIMMOMATIC_TRAILINGWINDOW
                        Trimmomatic window processing arguments, default:
                        TRAILING:3
  --trimmomatic_slidingwindow TRIMMOMATIC_SLIDINGWINDOW
                        Trimmomatic window processing arguments, default:
                        SLIDINGWINDOW:4:15
  --trimmomatic_quality TRIMMOMATIC_QUALITY
                        Trimmomatic quality encoding -phred33 or phred64
```

Example usage:
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

The specified assembler can be made through the `--method` option.
The full set of options are below.

```
usage: AAFTF assemble [-h] [-q] [--method METHOD] -o OUT [-w WORKDIR]
                      [-c cpus] [-m MEMORY] [-l LEFT] [-r RIGHT] [-v]
                      [--tmpdir TMPDIR] [--assembler_args ASSEMBLER_ARGS]
                      [--haplocontigs] [--pipe]

Run assembler on cleaned reads

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           Do not output warnings to stderr
  --method METHOD       Assembly method: spades, dipspades, megahit
  -o OUT, --out OUT     Output assembly FASTA
  -w WORKDIR, --workdir WORKDIR
                        assembly output directory
  -c cpus, --cpus cpus  Number of CPUs/threads to use.
  -m MEMORY, --memory MEMORY
                        Memory (in GB) setting for SPAdes. Default is 32
  -l LEFT, --left LEFT  Left (Forward) reads
  -r RIGHT, --right RIGHT
                        Right (Reverse) reads
  -v, --debug           Print Spades stdout to terminal
  --tmpdir TMPDIR       Assembler temporary dir
  --assembler_args ASSEMBLER_ARGS
                        Additional SPAdes/Megahit arguments
  --haplocontigs        For dipSPAdes take the haplocontigs file
  --pipe                AAFTF is running in pipeline mode
```

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
