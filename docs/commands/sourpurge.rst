=========
sourpurge
=========

**Aliases:** ``purge``

Purges contaminant contigs from an assembly using two independent signals: sourmash-based
taxonomic classification, and (if reads are provided) unusually low read coverage relative to the
bulk of the assembly.

Algorithm
=========

1. **Taxonomy screen.** Computes a per-contig sourmash sketch (``sourmash compute -k KMER
   --scaled=1000 --singleton``) and classifies each contig with
   ``sourmash lca classify`` against a reference LCA database. The database is either
   user-supplied (``--sourdb``) or auto-selected/downloaded via ``--sourdb_type``: ``gbk``
   (GenBank k=31, the default), ``gtdb``, or ``gtdbrep``. Any contig whose classified taxonomy
   does **not** include one of the phyla listed in ``-p/--phylum`` is marked for removal --
   contigs with *no* classification (``nomatch``) are kept (not penalized, since sourmash can't
   classify novel/divergent sequence).
2. **Low-coverage screen** (only run if ``-l/--left`` reads are supplied). Maps reads to the
   taxonomy-filtered assembly with ``bwa mem`` + ``samtools sort``, computes per-contig read
   coverage with ``samtools bedcov``, then computes the mean coverage of contigs at/above the
   assembly's **N50** contig length as a "true genome coverage" baseline. Any contig with
   coverage <= ``mincovpct%`` of that N50-contig average coverage is marked for removal (default
   5% -- i.e. contigs covered at less than 1/20th of the bulk-genome sequencing depth, typically
   evidence of a low-level contaminant or index-hopping/cross-contamination artifact rather than
   the target organism).
3. Contigs flagged by either screen are dropped; everything else is written to ``-o/--outfile``.
   A copy of the raw sourmash classification CSV is saved as
   ``{input_basename}.sourmash-taxonomy.csv`` next to the input file for review.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-p/--phylum``
     - *(required)*
     - Phylum/phyla to **keep** (e.g. ``Ascomycota``); anything classified outside this set is
       purged
   * - ``-k/--kmer``
     - 31
     - sourmash k-mer size (must match the database's k-mer size)
   * - ``--sourdb_type``
     - gbk
     - gbk (GenBank) / gtdb / gtdbrep sourmash LCA database
   * - ``-mc/--mincovpct``
     - 5
     - Minimum coverage, as a percent of the N50-contig average coverage, to avoid low-coverage
       removal

Invocation
==========

.. code-block:: text

    AAFTF sourpurge -i INPUT -o OUTFILE -p PHYLUM [PHYLUM ...]
                    [-l LEFT] [-r RIGHT] [-k KMER] [--sourdb PATH] [--sourdb_type {gbk,gtdb,gtdbrep}]
                    [-mc MINCOVPCT] [-c CPUS] [--AAFTF_DB DIR] [-w WORKDIR]
                    [--just-show-taxonomy] [-v] [--pipe]

``-i/--input``, ``-o/--outfile``, and ``-p/--phylum`` are required. Providing ``-l/--left`` (and
``-r/--right`` for paired data) enables the low-coverage screen; without reads, only the taxonomy
screen runs. ``--just-show-taxonomy`` prints classifications and exits without purging anything
(useful to sanity-check what phyla are present before choosing ``--phylum``).

Example
=======

.. code-block:: bash

    AAFTF sourpurge -c 24 -i genomes/STRAINX.vecscreen.fasta -o genomes/STRAINX.sourpurge.fasta \
        --left reads_filtered_1.fastq.gz --right reads_filtered_2.fastq.gz \
        --phylum Ascomycota --AAFTF_DB "$AAFTF_DB"

    # Preview classifications before deciding which phylum to keep
    AAFTF sourpurge -i genomes/STRAINX.vecscreen.fasta -o /dev/null \
        --phylum Ascomycota --just-show-taxonomy --AAFTF_DB "$AAFTF_DB"

Next step: :doc:`rmdup`.
