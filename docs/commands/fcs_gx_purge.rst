============
fcs_gx_purge
============

**Aliases:** ``ncbi_fcs-gx``, ``ncbi_fcs_gx``, ``gx``

Purges contaminant contigs from an assembly using NCBI's **FCS-GX** tool -- a genome
cross-contamination screen that classifies every contig/region against a comprehensive reference
database of taxonomy-labeled genomes, rather than the k-mer/sketch-based approach used by
:doc:`sourpurge`. This is an alternative to ``sourpurge`` for the taxonomic-contamination-purge
step of the pipeline; pick one or the other (or run both and compare).

Algorithm
=========

1. Runs NCBI's ``run_gx.py --fasta INPUT --tax-id TAXID --gx-db DB --out-dir WORKDIR``, which
   aligns/classifies assembly sequences against the FCS-GX reference database and writes a report
   named ``{input_basename}.{taxid}.fcs_gx_report.txt``.
2. Parses that report (tab-separated: ``seq_id start_pos end_pos seq_len action div
   agg_cont_cov top_tax_name``) and collects every sequence ID FCS-GX flagged for removal.
3. Writes all contigs **not** flagged to the output FASTA; a copy of the FCS-GX report is saved
   alongside the output as ``{outfile_basename}.fcs_gx-taxonomy.tsv`` for later review.

Requires the FCS-GX database
================================

Unlike ``sourpurge``, FCS-GX needs a very large (hundreds of GB), memory-mapped reference
database set up separately per NCBI's instructions
(https://github.com/ncbi/fcs/wiki/FCS-GX) -- AAFTF does not download or manage this database for
you. ``-d/--db`` must point at that database (default placeholder:
``/my_tmpfs/gxdb/all``, i.e. a fast local/tmpfs mount is strongly recommended given the database
size). ``fcs_gx_purge`` exits immediately with an error if ``{db}.gxi`` is not found.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-t/--taxid``
     - 4890
     - NCBI Taxonomy ID for the query organism (4890 = Ascomycota); FCS-GX uses this to decide
       what counts as "self" vs. contaminant
   * - ``-d/--db``
     - ``/my_tmpfs/gxdb/all``
     - Path to the pre-built FCS-GX database (must exist; large memory/SSD recommended)
   * - Contig removal criterion
     - any row present in the FCS-GX report
     - FCS-GX itself determines the ``action``/threshold logic internally; AAFTF removes every
       flagged sequence ID as-is

Invocation
==========

.. code-block:: text

    AAFTF fcs_gx_purge -i INPUT -o OUTFILE -d GXDB_PATH [-t TAXID] [-c CPUS]
                       [-w WORKDIR] [-v] [--pipe]

``-i/--input`` and ``-o/--outfile`` are required; ``-d/--db`` is required in practice (the run
aborts without a valid ``{db}.gxi``).

Example
=======

.. code-block:: bash

    AAFTF fcs_gx_purge -i genomes/STRAINX.vecscreen.fasta -o genomes/STRAINX.fcsgx.fasta \
        -d /fast_local/gxdb/all -t 4890 -c 32

Next step: :doc:`rmdup` (same as after :doc:`sourpurge`).
