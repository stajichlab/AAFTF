========
pipeline
========

No aliases.

Runs the entire raw-reads-to-clean-assembly workflow in a single command: trim -> mito (if paired
reads) -> filter -> assemble -> vecscreen -> sourpurge -> rmdup -> polish -> sort -> assess. See
:doc:`../workflow` for the full diagram and per-step input/output description.

Algorithm
=========

``pipeline`` is an orchestrator, not a separate implementation -- for each step it builds an
``argparse.Namespace`` from a filtered subset of the pipeline's own arguments (plus step-specific
required values, always including ``pipe=True`` so intermediate steps don't print "next command"
hints) and calls that submodule's ``run()`` function directly, in-process.

Before each step, ``pipeline`` checks whether that step's expected output file already exists; if
so, the step is **skipped** (with a log message) rather than re-run. This makes an interrupted
pipeline run resumable: simply re-invoke the same ``AAFTF pipeline`` command and it will pick up
after the last completed step. After each step, output existence is re-checked and the pipeline
aborts with an error if the expected output was not produced.

Mitochondrial assembly (:doc:`mito`) only runs when paired-end reads (``-r/--right``) are
provided; the resulting ``{basename}.mito.fasta``, if produced, is automatically passed to
:doc:`filter` as an extra ``--screen_local`` reference so mitochondrial reads are excluded from
the nuclear-genome read set rather than treated as contamination.

Steps intentionally **not** included in ``pipeline`` -- run these manually if needed:
:doc:`fcs_screen`, :doc:`fcs_gx_purge` (alternatives/complements to :doc:`vecscreen` /
:doc:`sourpurge`), :doc:`depth` (coverage QC of the final assembly), :doc:`fix_tbl` (post-FCS
annotation coordinate fixups), and :doc:`download` (run once, ahead of time, to populate
``$AAFTF_DB``).

Cutoffs / defaults
===================

``pipeline`` reuses each step's own defaults (:doc:`trim`, :doc:`mito`, :doc:`filter`,
:doc:`assemble`, :doc:`vecscreen`, :doc:`sourpurge`, :doc:`rmdup`, :doc:`polish`, :doc:`sort`,
:doc:`assess`) except where noted:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-m/--memory``
     - auto (75% of detected system RAM) if not set
     - Passed through to assemble/polish
   * - ``--method``
     - spades
     - Assembler method (spades / dipspades / megahit)
   * - ``-it/--iterations``
     - 5
     - Pilon polishing iterations
   * - ``-mc/--mincontiglen``
     - 500
     - Minimum contig length kept by rmdup and by the final sort step
   * - ``--mincovpct``
     - 5
     - sourpurge low-coverage removal threshold (percent of N50-contig coverage)

Invocation
==========

.. code-block:: text

    AAFTF pipeline -l LEFT [-r RIGHT] -o BASENAME -p PHYLUM [PHYLUM ...]
                   [-c CPUS] [-m MEMORY] [-ml MINLEN] [-it ITERATIONS]
                   [-mc MINCONTIGLEN] [--method {spades,dipspades,megahit}]
                   [-a ACCESSIONS ...] [-u URLS ...] [--sourdb PATH]
                   [--mincovpct PCT] [--AAFTF_DB DIR] [-w WORKDIR]
                   [--assembler_args ARG ...] [--tmpdir DIR] [-v] [--pipe]

``-l/--left``, ``-o/--out`` (basename), and ``-p/--phylum`` are required.

Example
=======

.. code-block:: bash

    export AAFTF_DB=~/lib/AAFTF_DB

    AAFTF pipeline \
        -l reads/STRAINX_R1.fq.gz -r reads/STRAINX_R2.fq.gz \
        -o STRAINX -c 24 -m 96 -it 5 \
        --phylum Ascomycota \
        --AAFTF_DB "$AAFTF_DB"

This produces, in sequence: ``STRAINX_1P.fastq.gz``/``STRAINX_2P.fastq.gz`` (trim),
``STRAINX.mito.fasta`` (mito, if paired), ``STRAINX_filtered_1.fastq.gz``/``_2.fastq.gz``
(filter), ``STRAINX.spades.fasta`` (assemble), ``STRAINX.vecscreen.fasta`` (vecscreen),
``STRAINX.sourpurge.fasta`` (sourpurge), ``STRAINX.rmdup.fasta`` (rmdup),
``STRAINX.polish.fasta`` (polish), and finally ``STRAINX.final.fasta`` with printed/``assess``
summary statistics.
