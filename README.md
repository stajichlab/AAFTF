### AAFTF - Automatic Assembly For The Fungi
*Jason Stajich and Jon Palmer*

Requirements
===================
- SPAdes - http://cab.spbu.ru/software/spades/
- Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic
- bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
- bwa - https://github.com/lh3/bwa
- Pilon - https://github.com/broadinstitute/pilon/wiki
- sourmash - https://sourmash.readthedocs.io/ (install via conda/pip)
- NCBI BLAST+ - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
- minimap2 - https://github.com/lh3/minimap2

Contributors
============
@hyphaltip
@nextgenusfs

Notes
===========
This is partially a python re-write of JAAWS - https://github.com/nextgenusfs/jaaws which was a unix shell 


Steps / Procedures
==================
1. trim - Read Trimmimg (Trimmomatic)
2. purge - Contaminant read purging - Bowtie2 matches to common contaminant (or user defined) datasets
3. assemble - Assembly (SPAdes)
4. vecscreen - Vector cleaning (BLAST vectorscreening against UniVec database)
5. sourpurge - SourMash removal of contaminanted contigs/scaffolds / (an alternative step is to use blobtools but is much slower)
7. rmdup - Funannotate clean for duplicated contig removal
8. pilon - Pilon for read-based assembly polishing
9. sort - sort and rename scaffolds by size
10. busco - BUSCO runs for summary statistics of gene content (code in progress)
