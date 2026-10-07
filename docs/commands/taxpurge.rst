========
taxpurge
========

Removes contigs that sourmash classifies outside the taxa you keep. Run it before :doc:`covpurge`,
so contaminant contigs do not skew covpurge's coverage baseline. It needs to know which phylum (or
other taxon) the organism belongs to; if you don't know, skip it -- the pipelines run it only when
given ``-p/--phylum``.

Algorithm
=========

1. Computes a per-contig sourmash sketch (``sourmash compute -k KMER --scaled=1000 --singleton``)
   and classifies each contig with ``sourmash lca classify`` against a reference LCA database. The
   database is either user-supplied (``--sourdb``) or the one installed by ``AAFTF database`` for
   ``--sourdb_type``: ``gbk`` -> ``sm_gbk`` (GenBank k=31, the default), ``gtdbrep`` ->
   ``sm_gtdbrep``, ``gtdb`` -> ``sm_gtdb`` (exit code 2 if it is not installed).
2. Only rows with status exactly ``found`` are classified: a contig is removed if none of the
   names given to ``-p/--phylum`` matches (exactly, case-sensitive) any rank of its lineage.
   Contigs with *no* classification (``nomatch``) are kept, since sourmash can't classify
   novel/divergent sequence.
3. The remaining contigs are written to ``-o/--outfile``, and a copy of the raw sourmash
   classification CSV is saved as ``{input_basename}.sourmash-taxonomy.csv`` next to the input.

Intermediate files and ``taxpurge.log`` are written to the work directory, which is kept when
given with ``-w``, with ``-v``, or when the run fails. Unless ``-q`` is given, the suggested next
command (:doc:`covpurge`) is printed.

Cutoffs / defaults
==================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-p/--phylum``
     - *(required)*
     - Phylum/phyla to **keep** (e.g. ``Ascomycota``); anything classified outside this set is
       removed
   * - ``-k/--kmer``
     - 31
     - sourmash k-mer size (must match the database's k-mer size)
   * - ``--sourdb_type``
     - gbk
     - gbk (GenBank) / gtdbrep / gtdb sourmash LCA database from ``AAFTF database``
       (ignored when ``--sourdb`` is given)

Invocation
==========

.. code-block:: text

    AAFTF taxpurge -i FASTA -o FASTA -p PHYLUM [PHYLUM ...] [-w DIR]
                   [--sourdb FILE] [-k KMER] [--sourdb_type {gbk,gtdbrep,gtdb}]
                   [--just-show-taxonomy] [-q] [-v]

``--just-show-taxonomy`` logs the distinct lineages found, saves
``{input_basename}.sourmash-taxonomy.csv`` next to the input, and returns without purging or writing
an output assembly (useful to check which phyla are present before choosing ``--phylum``).

Example
=======

.. code-block:: bash

    AAFTF taxpurge -i genomes/STRAINX.vecscreen.fasta -o genomes/STRAINX.taxpurge.fasta \
        --phylum Ascomycota

    # Preview classifications before deciding which phylum to keep
    AAFTF taxpurge -i genomes/STRAINX.vecscreen.fasta -o /dev/null \
        --phylum Ascomycota --just-show-taxonomy

Next step: :doc:`covpurge` (or :doc:`rmdup` without Illumina reads).
