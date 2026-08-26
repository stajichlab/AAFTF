======
filter
======

**Aliases:** ``filter_reads``, ``read_filter``

Removes reads that match known contaminant sequences (PhiX spike-in, UniVec, and any additional
sequences you specify) from trimmed FASTQ reads, before assembly.

Algorithm
=========

1. Builds (or reuses a cached) combined contamination FASTA (``contamdb.fa``) from:

   * the PhiX genome and UniVec (always included, downloaded from NCBI into ``$AAFTF_DB`` or the
     working directory on first use -- or pre-populated via :doc:`download`);
   * any ``-a/--screen_accessions`` GenBank accessions (fetched via NCBI eutils);
   * any ``-u/--screen_urls`` remote FASTA URLs;
   * any ``-s/--screen_local`` local FASTA files (e.g. a ``mito.fasta`` from :doc:`mito`, so
     mitochondrial reads aren't discarded as nuclear-genome contaminants).

2. Maps/kmer-matches the input reads against that combined database with one of four
   ``--aligner`` backends, then writes out only the **unmapped** (non-contaminant) reads:

   * **bbduk** (default) -- kmer-based filtering with ``k=27 hdist=1``, additionally matching
     built-in BBTools ``phix``/``artifacts``/``lambda`` reference sets. Paired-mode input is
     interleaved first (same PairStreamer workaround as :doc:`trim`) then de-interleaved after
     filtering.
   * **bowtie2** -- ``bowtie2 --very-sensitive``, builds a bowtie2 index of the contamination DB
     if missing/stale, output piped through ``samtools sort``.
   * **bwa** -- ``bwa mem`` against a bwa index of the contamination DB.
   * **minimap2** -- ``minimap2 -ax sr`` (short-read preset).

   For bowtie2/bwa/minimap2, reads are kept if flagged unmapped in the resulting sorted BAM
   (``samtools fastq -f 12`` for pairs -- both mates unmapped; ``-f 4`` for single-end).

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``--aligner``
     - bbduk
     - bbduk / bowtie2 / bwa / minimap2
   * - bbduk kmer size
     - 27
     - kmer length for contamination matching (``k=27``)
   * - bbduk hamming distance
     - 1
     - Allowed mismatches per kmer match (``hdist=1``)
   * - ``-m/--memory``
     - auto (60% of detected system RAM)
     - Max heap for bbduk (``-Xmx``)

Invocation
==========

.. code-block:: text

    AAFTF filter -l LEFT [-r RIGHT] [-o BASENAME] [-c CPUS]
                 [--aligner {bbduk,bowtie2,bwa,minimap2}] [-m MEMORY]
                 [-a ACCESSIONS ...] [-u URLS ...] [-s LOCAL_FASTA ...]
                 [--AAFTF_DB DIR] [-w WORKDIR] [-v] [--pipe]

**Output:** paired mode writes ``{basename}_filtered_1.fastq.gz`` /
``{basename}_filtered_2.fastq.gz``; single-end mode writes ``{basename}_filtered_U.fastq.gz``
(bbduk) or ``{basename}_filtered.fastq.gz`` (other aligners).

Example
=======

.. code-block:: bash

    AAFTF filter -c 16 --memory 64 --aligner bbduk \
        -o reads_trimmed/STRAINX \
        --left reads_trimmed/STRAINX_1P.fastq.gz --right reads_trimmed/STRAINX_2P.fastq.gz

    # also screen out a specific GenBank contaminant genome
    AAFTF filter -c 16 --aligner bbduk -o reads_trimmed/STRAINX \
        --left STRAINX_1P.fastq.gz --right STRAINX_2P.fastq.gz \
        -a NC_001422.1

Next step: :doc:`assemble`.
