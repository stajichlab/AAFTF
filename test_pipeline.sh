#!/bin/bash

module load ncbi-blast/2.7.1+
module load trimmomatic
module load java
module load bowtie2
module load python/3

if [ ! -d test ]; then
	mkdir test
	pushd test
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR522/005/SRR5223785/SRR5223785_1.fastq.gz
	curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR522/005/SRR5223785/SRR5223785_2.fastq.gz
	popd
fi

./scripts/AAFTF trim --left test/SRR5223785_1.fastq.gz --right test/SRR5223785_2.fastq.gz  --prefix Rant -c 24 --trimmomatic  $TRIMMOMATIC
./scripts/AAFTF filter -i test --prefix Rant --bowtie2 --paired
