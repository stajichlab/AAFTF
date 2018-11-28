### AAFTF - Automatic Assembly For The Fungi
*Jason Stajich and Jon Palmer*

Requirements
===================
- SPAdes - http://cab.spbu.ru/software/spades/
- BBTools - https://jgi.doe.gov/data-and-tools/bbtools/
- Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic (Optional)
- bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml (Optional)
- bwa - https://github.com/lh3/bwa
- Pilon - https://github.com/broadinstitute/pilon/wiki
- sourmash - https://sourmash.readthedocs.io/ (install via conda/pip)
- NCBI BLAST+ - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
- minimap2 - https://github.com/lh3/minimap2

Authors
============
* Jason Stajich [@hyphaltip](https://github.com/hyphaltip) - http://lab.stajich.org
* Jon Palmer [@nextgenusfs](https://github.com/nextgenusfs) - https://twitter.com/jonpalmer2013

Notes
===========
This is partially a python re-write of [JAAWS](https://github.com/nextgenusfs/jaaws) which was a unix shell based cleanup and assembly tool written by Jon.


Steps / Procedures
==================
1. trim                Trim FASTQ input reads
2. filter              Filter contaminanting reads
3. assemble            Assemble reads
4. vecscreen           Vector and Contaminant Screening of assembled contigs
5. sourpurge           Purge contigs based on sourmash results
6. rmdup               Remove duplicate contigs
7. pilon               Polish contig sequences with Pilon
8. sort                Sort contigs by length and rename FASTA headers
9. assess              Assess completeness of genome assembly
10. pipeline            Run AAFTF pipelinetrim - Read Trimmimg (Trimmomatic)
