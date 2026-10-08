========
taxpurge
========

Removes contigs that sourmash classifies outside the taxa you keep. Run it before :doc:`covpurge`,
so contaminant contigs do not skew covpurge's coverage baseline. It needs to know which phylum (or
other taxon) the organism belongs to; if you don't know, skip it -- the pipelines run it only when
given ``-p/--phylum``.

Algorithm
=========

1. Sketches each contig (``sourmash sketch dna -p k=KMER,scaled=1000 --singleton``).
2. Searches the contig sketches against each ``--sourdb`` database in turn, at that database's
   resolution (``scaled``), with ``sourmash scripts fastmultigather`` from the `branchwater plugin
   <https://github.com/sourmash-bio/sourmash_plugin_branchwater>`_. A contig matches a database
   genome when they share at least ``--min_hashes`` sketch hashes, i.e. about ``--min_hashes`` x
   ``scaled`` bp of sequence.
3. Classifies the matches of each database with ``sourmash tax genome``. A contig matched in
   several databases takes the classification with the most matched sequence, so a contig of the
   target organism is not called a contaminant because a (contaminated) reference genome contains
   a little of it.
4. Removes classified contigs whose lineage contains none of the ``-p/--phylum`` names (exact,
   case-sensitive match at any rank; GTDB rank prefixes such as ``p__`` are optional). Unclassified
   contigs are kept, since sourmash can't classify novel or divergent sequence.
5. Writes the remaining contigs to ``-o/--outfile`` and saves the chosen classifications as
   ``{input_basename}.sourmash-taxonomy.csv`` next to the input; contigs missing from it were not
   classified.

Databases
=========

Install them with ``AAFTF database NAME`` (which also downloads the lineage CSV), or ``AAFTF
database taxpurge`` for the default set. All are k=31.

.. list-table::
   :header-rows: 1
   :widths: 18 52 15 15

   * - ``--sourdb`` name
     - Contents
     - Size
     - ``scaled``
   * - ``gbk-fungi``
     - GenBank 2025.04 fungal genomes (21,135)
     - 5.3 GB
     - 1000
   * - ``gbk-fungi-lite`` (default)
     - GenBank 2025.01 eukaryote collection, fungi (~6,300 genomes)
     - 170 MB
     - 10000
   * - ``gtdb-reps``
     - GTDB rs226 representative genomes: bacteria and archaea (143,384)
     - 3.9 GB
     - 1000
   * - ``gtdb-reps-lite`` (default)
     - the same genomes at lower resolution
     - 447 MB
     - 10000
   * - ``human-hg38`` (default)
     - human reference genome hg38 (taxpurge writes its one-line lineage CSV itself)
     - 55 MB
     - 1000

Pick one fungal and one GTDB database (any combination), and add ``human-hg38`` to catch human
contamination. Search the target organism's group too (the fungi databases for a fungal
assembly): with bacteria alone, contigs of your organism can be called bacterial, because some
bacterial genome assemblies contain sequence from other organisms. ``--sourdb`` also takes paths
of other sourmash sketch collections (their ``scaled`` is read from the collection) or branchwater
RocksDB indexes (assumed ``scaled=1000``), with their lineage CSVs in ``--taxonomy_csv``.

**Lite or full.** At ``scaled=10000`` a contig needs about ``--min_hashes`` x 10 kb of sequence
matching one genome to be classified, at ``scaled=1000`` about ``--min_hashes`` x 1 kb. In a test
(yeast assembly with *E. coli*, *Cutibacterium acnes* and human chr21 contigs mixed in, 4 threads):

.. list-table::
   :header-rows: 1
   :widths: 40 18 22 20

   * - Databases
     - Time / memory
     - Contaminants removed
     - Yeast wrongly removed
   * - ``gbk-fungi-lite gtdb-reps-lite human-hg38`` (default)
     - 43 s / 1.1 GB
     - bacteria 38/315 (none under 5 kb), human 33/40
     - 0/698
   * - ``gbk-fungi gtdb-reps`` (bacteria only, measured with ``--min_hashes 3``)
     - 3.7 min / 7.5 GB [#]_
     - bacteria 143/315 (about 1 in 7 under 5 kb)
     - 2/698 (3-4 kb)

.. [#] Peak with ``gbk-fungi`` searched from its RocksDB index (see below); with its sketch file
   (``AAFTF database --no-build``), memory is higher.

Use the full databases when short contaminant contigs matter and memory allows; contaminants on
short contigs are otherwise better caught by :doc:`covpurge` or :doc:`fcs_gx_purge`.

**Memory and RocksDB indexes.** Each database is searched separately, so memory is that of the
largest one. ``AAFTF database`` installs each database as a branchwater RocksDB index, which is
searched from disk and uses much less memory: ``gtdb-reps`` and ``gtdb-reps-lite`` are downloaded as
published indexes; ``gbk-fungi``, ``gbk-fungi-lite`` and ``human-hg38`` are downloaded as sketch
files, indexed locally (``gbk-fungi``: about 30 min, 6 GB of memory, 8.6 GB of disk) and the
sketch files removed. ``AAFTF database --no-build`` keeps the sketch files of those three
instead; ``taxpurge`` then searches the sketch file, loading it into memory. With all three
default databases as indexes, the test above took 37 s and 1.2 GB; with ``gbk-fungi`` from its
index plus ``gtdb-reps-lite`` and ``human-hg38``, 39 s and 1.1 GB. Other RocksDB indexes can be
passed to ``--sourdb`` as paths, with their lineage CSVs in ``--taxonomy_csv``.

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
   * - ``--sourdb``
     - gbk-fungi-lite gtdb-reps-lite human-hg38
     - Databases to search (see above), or paths of sourmash sketch collections / RocksDB indexes
   * - ``--taxonomy_csv``
     - none
     - Lineage CSVs for databases given as paths
   * - ``--min_hashes``
     - 2
     - Sketch hashes a contig must share with a database genome to match it (x ``scaled`` = bp)
   * - ``-k/--kmer``
     - 31
     - sourmash k-mer size (must be in the databases)
   * - ``-c/--cpus``
     - 1
     - Threads for ``fastmultigather``

Invocation
==========

.. code-block:: text

    AAFTF taxpurge -i FASTA -o FASTA -p PHYLUM [PHYLUM ...] [-w DIR]
                   [--sourdb DB [DB ...]] [--taxonomy_csv CSV [CSV ...]] [-k KMER]
                   [--min_hashes INT] [-c INT] [--just-show-taxonomy] [-q] [-v]

``--just-show-taxonomy`` logs how many contigs fall in each phylum, saves
``{input_basename}.sourmash-taxonomy.csv`` next to the input, and returns without purging or writing
an output assembly (useful to check which phyla are present before choosing ``--phylum``).

Example
=======

.. code-block:: bash

    AAFTF database taxpurge    # the default databases, once (~0.7 GB)
    AAFTF taxpurge -c 8 -i genomes/STRAINX.vecscreen.fasta -o genomes/STRAINX.taxpurge.fasta \
        --phylum Ascomycota

    # Full-resolution databases: shorter contaminant contigs caught, more memory
    AAFTF database gbk-fungi gtdb-reps human-hg38
    AAFTF taxpurge -c 8 -i genomes/STRAINX.vecscreen.fasta -o genomes/STRAINX.taxpurge.fasta \
        --phylum Ascomycota --sourdb gbk-fungi gtdb-reps human-hg38

    # Preview classifications before deciding which phylum to keep
    AAFTF taxpurge -i genomes/STRAINX.vecscreen.fasta -o /dev/null \
        --phylum Ascomycota --just-show-taxonomy

Next step: :doc:`covpurge` (or :doc:`rmdup` without Illumina reads).
