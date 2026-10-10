============
fcs_gx_purge
============

Purges contaminant contigs from an assembly using NCBI's **FCS-GX** tool -- a genome
cross-contamination screen that classifies every contig/region against a comprehensive reference
database of taxonomy-labeled genomes, rather than the k-mer/sketch-based approach used by
:doc:`taxpurge`. This is an alternative to ``taxpurge`` for the taxonomic-contamination-purge
step of the pipeline; pick one or the other (or run both and compare).

Algorithm
=========

1. Runs NCBI's ``run_gx.py --fasta INPUT --tax-id TAXID --gx-db DB --out-dir WORKDIR``, which
   aligns/classifies assembly sequences against the FCS-GX reference database and writes a report
   named ``{input_basename}.{taxid}.fcs_gx_report.txt``.
2. Parses that report (tab-separated: ``seq_id start_pos end_pos seq_len action div
   agg_cont_cov top_tax_name``). Only rows whose action is ``EXCLUDE`` mark a sequence for
   removal; ``TRIM``, ``FIX``, ``REVIEW`` and other actions are logged as warnings and the
   sequence is **kept** (partial-contig edits are not applied -- review them manually).
3. Writes all contigs not excluded to the output FASTA; a copy of the FCS-GX report is saved
   alongside the output as ``{outfile_basename}.fcs_gx-taxonomy.tsv`` for later review.

``run_gx.py`` stderr goes to ``fcs_gx.log`` in the work directory, which also holds
``fcs_gx_purge.log`` and is kept when given with ``-w``, with ``-v``, or when the run fails. If
``run_gx.py`` exits non-zero or writes no report, the command fails with a ``RuntimeError``
(exit code 1). Unless ``-q`` is given, the suggested next command (:doc:`rmdup`) is printed.

Requires the FCS-GX database
================================

Unlike ``taxpurge``, FCS-GX needs a very large (hundreds of GB), memory-mapped reference
database set up separately per NCBI's instructions
(https://github.com/ncbi/fcs/wiki/FCS-GX) -- AAFTF does not download or manage this database for
you. ``-d/--db`` must point at that database (default placeholder:
``/my_tmpfs/gxdb/all``, i.e. a fast local/tmpfs mount is strongly recommended given the database
size). ``fcs_gx_purge`` stops immediately (``FileNotFoundError``, exit code 2) if ``{db}.gxi``
is not found.

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
     - ``EXCLUDE`` action
     - FCS-GX decides the ``action`` internally; AAFTF drops only ``EXCLUDE`` sequences and
       warns about (but keeps) everything else

Invocation
==========

.. code-block:: text

    AAFTF fcs_gx_purge -i FASTA -o FASTA [-w DIR] [-t TAXID] [-d DB] [-q] [-v]

``-i/--input`` and ``-o/--outfile`` are required; ``-d/--db`` must point at a real database
unless it lives at the default path (the run aborts without a valid ``{db}.gxi``).

Example
=======

.. code-block:: bash

    AAFTF fcs_gx_purge -i genomes/STRAINX.vecscreen.fasta -o genomes/STRAINX.fcsgx.fasta \
        -d /fast_local/gxdb/all -t 4890

Next step: :doc:`covpurge` (same as after :doc:`taxpurge`).
