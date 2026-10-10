====
sort
====

Sorts final assembly contigs by length (longest to shortest) and renames FASTA headers to a
clean, consistent scheme -- typically the last cosmetic step before assessment/submission.

Algorithm
=========

Reads all sequences, drops any shorter than ``-ml/--minlen``, de-duplicates by header (first
occurrence wins), sorts the remainder by sequence length descending, and writes them out renamed
``{name}_1``, ``{name}_2``, ... in that order (original header text is discarded). No work
directory is used; ``./sort.log`` is written only with ``-v``.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-ml/--minlen``
     - 0 (no filtering)
     - Minimum contig length retained
   * - ``-n/--name``
     - ``scaffold``
     - Basename prefix for renamed FASTA headers (``scaffold_1``, ``scaffold_2``, ...)

Invocation
==========

.. code-block:: text

    AAFTF sort -i FASTA -o FASTA [-ml BP] [-n NAME] [-q] [-v]

``-i/--input`` (alias ``--infile``) and ``-o/--out`` (alias ``--output``) are required;
``-n/--name`` also accepts ``--basename``.

Example
=======

.. code-block:: bash

    AAFTF sort -i genomes/STRAINX.rmdup.fasta -o genomes/STRAINX.final.fasta -n STRAINX

Next step: :doc:`assess` (and, optionally, :doc:`depth`).
