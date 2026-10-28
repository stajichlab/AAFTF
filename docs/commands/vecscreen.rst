=========
vecscreen
=========

Trims and splits adaptor/vector sequence out of assembled contigs, applying the rules NCBI's
FCS-adaptor uses, with one BLASTN search against a database built for fungal assemblies. It is an
alternative to :doc:`fcs_screen`, which runs NCBI's FCS-adaptor itself in a container: run one or
the other, not both.

Database
========

``vecscreen`` searches the ``contamdb`` BLAST database, shared with :doc:`filter` and built by
``AAFTF database contamdb`` (see :doc:`database`) from NCBI's adaptors for eukaryote assemblies,
PhiX and the BBTools ``lambda.fa.gz`` and ``adapters.fa``. It holds no cloning vectors, whose
fungal-derived parts split genuine fungal contigs, so ``vecscreen`` removes sequencing adaptors and
spike-ins, not vector inserts; add a vector's sequence to ``contamdb.fa`` to screen for it. The run stops with exit code 2
if ``contamdb.fa`` is missing; its BLAST database is rebuilt first if it is missing or older than
``contamdb.fa`` (e.g. after you added sequences to it).

Algorithm
=========

The assembly is searched once with blastn using NCBI's VecScreen settings (``-reward 1 -penalty -5
-gapopen 3 -gapextend 3 -dust yes -evalue 700 -searchsp 1750000000000``: the BLAST library's
VecScreen preset, given explicitly because the ``blastn`` command line has no ``-task vecscreen``).
Each hit's raw score is then judged as FCS-adaptor does:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Hit
     - Minimum score
     - Action
   * - Terminal, reaching the first/last base
     - 16
     - Trimmed off (ACTION_TRIM)
   * - Terminal (within 50 bp of an end)
     - 19
     - Trimmed off, together with the bases between it and the end (ACTION_TRIM)
   * - Internal (more than 50 bp from both ends)
     - 30 (Strong)
     - Cut out, splitting the contig into ``split<N>_<contig>`` pieces

Weak and moderate internal matches (e.g. short k-mer matches in conserved genes) are ignored, so
they never fragment chromosomes. Overlapping hits, and hits less than 50 bp apart, are merged into
one region. A contig with more than 50% of its length in vector regions, and any contig or piece
shorter than 200 bp after trimming, is dropped (ACTION_EXCLUDE).

The cleaned assembly is written to ``-o/--outfile``. The BLAST report goes to the work directory,
together with ``vecscreen.log`` (the work directory is kept when given with ``-w``, with ``-v``, or
when the run fails).

Invocation
==========

.. code-block:: text

    AAFTF vecscreen -i FASTA -o FASTA [-w DIR] [-c INT] [-q] [-v]

``-i/--input`` (alias ``--infile``; assembly FASTA) and ``-o/--outfile`` are required. Unless
``-q`` is given, the suggested next command (:doc:`taxpurge`) is printed at the end.

Example
=======

.. code-block:: bash

    AAFTF database contamdb   # once
    AAFTF vecscreen -c 16 -i genomes/STRAINX.spades.fasta -o genomes/STRAINX.vecscreen.fasta

Next step: :doc:`taxpurge` (sourmash-based taxonomic purge). :doc:`fcs_gx_purge` (NCBI FCS-GX) is
an alternative run separately, on a large-memory node, on finished assemblies.
