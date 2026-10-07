=====
rmdup
=====

Removes contigs that are wholly or largely redundant with a larger contig elsewhere in the same
assembly -- a common artifact of assemblers producing overlapping/haplotype-divergent contigs for
the same genomic region.

Algorithm
=========

1. Computes assembly N50 and **N75** contig lengths.
2. Iterates contigs from shortest to longest. Each contig shorter than ``-ml/--minlen`` is
   dropped outright without alignment. Each remaining contig shorter than the N75 length (or, with
   ``--exhaustive``, every contig regardless of length) is aligned against all *longer* contigs
   in the assembly with ``minimap2 -x asm5 -N5``.
3. A contig is discarded as a duplicate if any alignment is strictly greater than **both**
   ``-pid/--percent_id`` sequence identity **and** ``-pcov/--percent_cov`` query coverage
   (identity = matched bases / alignment length; coverage = alignment length / query length).
4. All non-duplicate, non-too-short contigs are written to ``-o/--out``. Temporary query/reference
   FASTA files and ``rmdup.log`` go to the work directory (kept when given with ``-w``, with
   ``-v``, or when the run fails). Unless ``-q`` is given, the suggested next command
   (:doc:`sort`) is printed.

Limiting the search to contigs <= N75 by default is a performance optimization: it assumes the
handful of large, low-copy scaffolds that make up the bulk of assembly length are unlikely to be
near-exact duplicates of each other, so only the many small "extra" contigs actually need the
exhaustive per-contig minimap2 search. Use ``--exhaustive`` if you suspect large-scale
duplication (e.g. after a diploid/haplotig-aware assembly).

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-pid/--percent_id``
     - 95
     - Minimum minimap2 alignment identity (%) to call a contig a duplicate
   * - ``-pcov/--percent_cov``
     - 95
     - Minimum query coverage (%) of the alignment to call a contig a duplicate
   * - ``-ml/--minlen``
     - 500
     - Contigs shorter than this are dropped unconditionally (too short to be reliably useful)
   * - Search scope
     - contigs < N75
     - ``--exhaustive`` checks every contig regardless of length

Invocation
==========

.. code-block:: text

    AAFTF rmdup -i FASTA -o FASTA [-c INT] [-w DIR] [-pid PCT] [-pcov PCT]
                [-ml BP] [--exhaustive] [-q] [-v]

``-i/--input`` and ``-o/--out`` are required.

Example
=======

.. code-block:: bash

    AAFTF rmdup -c 16 -i genomes/STRAINX.covpurge.fasta -o genomes/STRAINX.rmdup.fasta

    # more thorough (checks every contig, not just those <= N75)
    AAFTF rmdup -c 16 --exhaustive -i genomes/STRAINX.covpurge.fasta -o genomes/STRAINX.rmdup.fasta

Next step: :doc:`sort` (optionally run :doc:`polish` first).
