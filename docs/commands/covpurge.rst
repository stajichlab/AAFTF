========
covpurge
========

Removes contigs whose Illumina read coverage is very low compared to the bulk of the assembly --
typically low-level contaminants or index-hopping/cross-contamination artifacts rather than the
target organism. It needs no taxonomy, so it can run when the organism's phylum is unknown. Run it
after :doc:`taxpurge` when you use both.

Algorithm
=========

1. Maps the reads to the assembly with ``bwa mem`` + ``samtools sort`` and computes per-contig
   read coverage with ``samtools bedcov``.
2. Computes the mean coverage of contigs at/above the assembly's **N50** contig length as the
   "true genome coverage" baseline.
3. Removes any contig with coverage <= ``mincovpct%`` of that baseline (default 5%, i.e. contigs
   covered at less than 1/20th of the bulk-genome sequencing depth). If no N50-length contig has
   coverage data, no contig is removed and a warning is logged.

The BAM of a previous run in the same work directory (``-w``) is reused. Intermediate files and
``covpurge.log`` are written to the work directory, which is kept when given with ``-w``, with
``-v``, or when the run fails. Unless ``-q`` is given, the suggested next command (:doc:`rmdup`) is
printed.

Cutoffs / defaults
==================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-1/--read1``
     - *(required)*
     - Illumina read 1 (or single-end) FASTQ; add ``-2/--read2`` for paired reads
   * - ``-mc/--mincovpct``
     - 5
     - Minimum coverage, as a percent of the N50-contig average coverage, to avoid removal

Invocation
==========

.. code-block:: text

    AAFTF covpurge -i FASTA -o FASTA -1 FASTQ [-2 FASTQ] [-w DIR] [-mc PCT]
                   [-c INT] [-q] [-v]

Example
=======

.. code-block:: bash

    AAFTF covpurge -c 24 -i genomes/STRAINX.taxpurge.fasta -o genomes/STRAINX.covpurge.fasta \
        -1 reads_filtered_1.fastq.gz -2 reads_filtered_2.fastq.gz

Next step: :doc:`rmdup`.
