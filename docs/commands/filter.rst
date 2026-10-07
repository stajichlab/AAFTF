======
filter
======

Removes reads that match known contaminant sequences (PhiX spike-in, UniVec, and any additional
sequences you specify) from trimmed FASTQ reads, before assembly.

Algorithm
=========

1. Builds (or reuses a cached) combined contamination FASTA (``contamdb.fa``) from:

   * the PhiX genome and UniVec (always included; they are never downloaded by ``filter`` -- install
     them first with ``AAFTF database phix univec``, see :doc:`database`; a missing database
     raises an error, exit code 2);
   * the Illumina sequencing-artifact and lambda phage (a common spike-in) FASTAs bundled with
     BBTools (``sequencing_artifacts.fa.gz`` and ``lambda.fa.gz`` in the ``resources/`` folder
     next to ``bbduk.sh``; always included, error if BBTools is not found);
   * any ``-a/--screen_accessions`` GenBank accessions (downloaded via NCBI eutils);
   * any ``-u/--screen_urls`` remote FASTA URLs (downloaded);
   * any ``-s/--screen_local`` local FASTA files (e.g. a ``mito.fasta`` from :doc:`mito`, so
     mitochondrial reads aren't discarded as nuclear-genome contaminants).

2. Maps/kmer-matches the input reads against that combined database with one of four
   ``--aligner`` backends, then writes out only the **unmapped** (non-contaminant) reads:

   * **bbduk** (default) -- kmer-based filtering with ``k=27 hdist=1``. Paired-mode input is
     interleaved first (same PairStreamer workaround as :doc:`trim`) then de-interleaved after
     filtering.
   * **bowtie2** -- ``bowtie2 --very-sensitive``, builds a bowtie2 index of the contamination DB
     if missing/stale, output piped through ``samtools sort``.
   * **bwa** -- ``bwa mem`` against a bwa index of the contamination DB.
   * **minimap2** -- ``minimap2 -ax sr`` (short-read preset).

   Any other ``--aligner`` value raises an error. Every backend screens the same sequences: the
   whole combined database above.

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
     - 8 (GB)
     - Max Java heap for bbduk (``-Xmx``)

Invocation
==========

.. code-block:: text

    AAFTF filter -1 FASTQ [-2 FASTQ] [-o PREFIX] [-c INT]
                 [--aligner {bbduk,bowtie2,bwa,minimap2}] [-m GB]
                 [-a ACCESSION ...] [-u URL ...] [-s FASTA ...]
                 [-w DIR] [-q] [-v]

If ``-o/--out`` is not given, the prefix is derived from ``--read1``'s file name. Without
``-w/--workdir`` a temporary directory is used and removed afterwards (kept with ``-v`` or on
failure); the step log is written to ``<workdir>/filter.log``. In a reused ``-w`` directory,
``contamdb.fa`` is rebuilt whenever its list of sources (or any source file) changes, and the
aligner index and alignment are redone when they are older than ``contamdb.fa``.

**Output:** paired mode writes ``{basename}_filtered_1.fastq.gz`` /
``{basename}_filtered_2.fastq.gz``; single-end mode writes ``{basename}_filtered_1.fastq.gz``
(with every aligner).

Example
=======

.. code-block:: bash

    AAFTF filter -c 16 --memory 64 --aligner bbduk \
        -o reads_trimmed/STRAINX \
        --read1 reads_trimmed/STRAINX_1P.fastq.gz --read2 reads_trimmed/STRAINX_2P.fastq.gz

    # also screen out a specific GenBank contaminant genome
    AAFTF filter -c 16 --aligner bbduk -o reads_trimmed/STRAINX \
        --read1 STRAINX_1P.fastq.gz --read2 STRAINX_2P.fastq.gz \
        -a NC_001422.1

Next step: :doc:`assemble`.
