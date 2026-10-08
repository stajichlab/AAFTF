=========
Pipelines
=========

Three subcommands run the whole raw-reads-to-clean-assembly workflow in a single command:

* ``pipeline_short`` -- Illumina reads: trim -> filter -> assemble -> vecscreen -> [taxpurge] ->
  covpurge -> rmdup -> sort -> assess.
* ``pipeline_long`` -- long reads (ONT/PacBio) only: assemble (flye) -> vecscreen ->
  [taxpurge] -> rmdup -> polish (racon) -> sort -> assess.
* ``pipeline_hybrid`` -- Illumina + long reads: trim -> filter -> assemble (flye or unicycler) ->
  vecscreen -> [taxpurge] -> covpurge -> rmdup -> polish -> sort -> assess.

:doc:`taxpurge` runs only when ``-p/--phylum`` is given, so the pipelines also work when the
organism's taxonomy is unknown (the contigs are then screened by :doc:`vecscreen`,
:doc:`covpurge` and :doc:`rmdup` only). :doc:`covpurge` needs Illumina reads, so ``pipeline_long``
does not run it. Polishing always runs after the cleaning steps, so no time is spent polishing
contigs that are then discarded.

See :doc:`../workflow` for the full diagram and per-step input/output description.

.. note::

   ``pipeline_short`` was called ``pipeline`` before; there is no alias for the old name.

pipeline_short
==============

Algorithm
---------

Each pipeline is an orchestrator, not a separate implementation -- it calls each submodule's
``run()`` function directly, in-process. Every step starts from exactly the defaults it has on its
own command line (read from the ``AAFTF <step>`` parser, so the pipeline cannot drift from them);
only the options given to the pipeline and the file names that chain the steps together override
them. The steps' "next command" hints are not shown, since the pipeline runs the next step itself.

Before each step, the pipeline checks whether that step's expected output file already exists; if
so, the step is **skipped** (with a log message) rather than re-run. This makes an interrupted
pipeline run resumable: simply re-invoke the same ``AAFTF pipeline_short`` command and it will pick up
after the last completed step. After each step, output existence is re-checked and the pipeline
aborts with an error if the expected output was not produced.

Steps intentionally **not** included in ``pipeline_short`` (marked "(Optional)" in ``AAFTF -h``) -- run
these manually if needed: :doc:`mito` (to screen mitochondrial reads out, run it before
:doc:`filter` and pass its output to ``filter --screen_local``), :doc:`fcs_screen`, :doc:`fcs_gx_purge` (alternatives/complements to :doc:`vecscreen` /
:doc:`taxpurge`), :doc:`polish` (recommended only with long reads; polishing a short-read-only
assembly with the same reads rarely helps), :doc:`depth` (coverage QC of the final assembly), :doc:`fix_tbl` (post-FCS
annotation coordinate fixups), and :doc:`database` (run once, ahead of time, to populate
``$AAFTF_DB``).

Cutoffs / defaults
------------------

``pipeline_short`` uses each step's own defaults (:doc:`trim`, :doc:`filter`,
:doc:`assemble`, :doc:`vecscreen`, :doc:`taxpurge`, :doc:`covpurge`, :doc:`rmdup`, :doc:`sort`,
:doc:`assess`) except where noted:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-m/--memory``
     - each step's own default
     - Passed to every step that has ``-m/--memory`` (trim, filter, assemble)
   * - ``-c/--cpus``, ``-w/--workdir``, ``-v/--verbose``, ``-q/--quiet``
     - each step's own default
     - Passed to every step that has the option; ``assemble`` gets its own subfolder
       ``<workdir>/assemble_<method>`` (an assembler treats an existing folder as its own earlier run)
   * - ``-a/--screen_accessions``, ``-u/--screen_urls``
     - none
     - Extra sequences screened out of the reads by filter
   * - ``--tmpdir``, ``--assembler_args``
     - none
     - Passed to assemble
   * - ``-p/--phylum``, ``--sourdb``
     - none (taxpurge not run); ``--sourdb``: ``gbk-fungi-lite gtdb-reps-lite human-hg38``
     - Phyla kept, and the sourmash databases (see :doc:`taxpurge`), for taxpurge
   * - ``--method``
     - spades
     - Assembler method (spades / megahit / unicycler)
   * - ``-ml/--minlen``
     - 75
     - Minimum read length kept by trim
   * - ``-mc/--mincontiglen``
     - 500
     - Minimum contig length kept by rmdup and by the final sort step
   * - ``--mincovpct``
     - 5
     - covpurge low-coverage removal threshold (percent of N50-contig coverage)

Invocation
----------

.. code-block:: text

    AAFTF pipeline_short -1 FASTQ [-2 FASTQ] -o BASENAME [-p PHYLUM [PHYLUM ...]]
                         [-c CPUS] [-m MEMORY] [-ml MINLEN]
                         [-mc MINCONTIGLEN] [--method {spades,megahit,unicycler}]
                         [-a ACCESSIONS ...] [-u URLS ...] [--sourdb DB ...]
                         [--mincovpct PCT] [-w WORKDIR]
                         [--assembler_args ARG ...] [--tmpdir DIR] [-q] [-v]

``-1/--read1`` and ``-o/--out`` (basename) are required.

Example
-------

.. code-block:: bash

    export AAFTF_DB=~/lib/AAFTF_DB

    AAFTF pipeline_short \
        -1 reads/STRAINX_R1.fq.gz -2 reads/STRAINX_R2.fq.gz \
        -o STRAINX -c 24 -m 96 \
        --phylum Ascomycota

This produces, in sequence: ``STRAINX_1P.fastq.gz``/``STRAINX_2P.fastq.gz`` (trim),
``STRAINX_filtered_1.fastq.gz``/``_2.fastq.gz``
(filter), ``STRAINX.spades.fasta`` (assemble), ``STRAINX.vecscreen.fasta`` (vecscreen),
``STRAINX.taxpurge.fasta`` (taxpurge), ``STRAINX.covpurge.fasta`` (covpurge),
``STRAINX.rmdup.fasta`` (rmdup), and finally
``STRAINX.final.fasta`` (sort) with printed/``assess``
summary statistics.

pipeline_long
=============

Assembles long reads with :doc:`assemble` ``--method flye``, runs :doc:`vecscreen` and
:doc:`rmdup`, polishes the cleaned contigs once with :doc:`polish` ``--method racon`` using the same
reads, then runs :doc:`sort` and :doc:`assess`. With ``-p/--phylum``, :doc:`taxpurge` also runs
between vecscreen and rmdup. :doc:`covpurge` is not run, since it needs Illumina reads for contig
coverage. Steps are run, defaulted and skipped (resumable) exactly as in ``pipeline_short``.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-lr/--longreads``, ``-o/--out`` (required)
     -
     - Long-read FASTQ and output file prefix
   * - ``--longread_type``
     - nano-hq
     - nano-raw / nano-hq / pacbio-raw / pacbio-hifi (flye)
   * - ``--genome_size``
     - estimated by flye
     - Estimated genome size for flye, e.g. ``40m``
   * - ``-p/--phylum``
     - none (taxpurge not run)
     - Phyla whose contigs taxpurge keeps
   * - ``--sourdb``
     - from ``AAFTF database``
     - sourmash databases for taxpurge (see :doc:`taxpurge`)
   * - ``--assembler_args``
     - none
     - Passed to assemble
   * - ``-m/--memory``
     - each step's own default
     - Passed to polish
   * - ``-mc/--mincontiglen``
     - 500
     - Minimum contig length kept by rmdup and sort

.. code-block:: text

    AAFTF pipeline_long -lr FASTQ -o PREFIX
                        [--longread_type {nano-raw,nano-hq,pacbio-raw,pacbio-hifi}]
                        [--genome_size SIZE] [-p PHYLUM [PHYLUM ...]] [--sourdb DB ...]
                        [-w DIR] [--assembler_args ARG]
                        [-mc BP] [-c INT] [-m GB] [-q] [-v]

.. code-block:: bash

    AAFTF pipeline_long -lr reads/STRAINX_ont.fq.gz -o STRAINX \
        --genome_size 40m -c 24

This produces ``STRAINX.flye.fasta`` (assemble), ``STRAINX.vecscreen.fasta``,
``STRAINX.rmdup.fasta``, ``STRAINX.racon.fasta`` (polish) and ``STRAINX.final.fasta``, then prints the
``assess`` statistics.

pipeline_hybrid
===============

Trims and filters the Illumina reads (:doc:`trim`, :doc:`filter`), assembles the long reads with
flye, runs :doc:`vecscreen`, :doc:`taxpurge` (only with ``-p/--phylum``), :doc:`covpurge` (with the
Illumina reads) and :doc:`rmdup`, polishes
the cleaned contigs with :doc:`polish` ``--method racon`` (long reads) and then with the filtered
Illumina reads -- ``polypolish`` for paired reads, ``pypolca`` for single-end -- and runs
:doc:`sort` and :doc:`assess`. With ``--method unicycler``, unicycler assembles the filtered Illumina reads together
with the long reads and both polish steps are skipped, since Unicycler polishes its own assembly.
Steps are run, defaulted and skipped (resumable) exactly as in ``pipeline_short``.

Options are those of ``pipeline_short`` (except that ``--method`` is ``flye`` (default) or
``unicycler``) plus ``-lr/--longreads`` (required), ``--longread_type`` and ``--genome_size`` as
in ``pipeline_long``. ``-m/--memory`` is passed to trim, filter, assemble and polish.

.. code-block:: text

    AAFTF pipeline_hybrid -1 FASTQ -lr FASTQ -o PREFIX [-p PHYLUM [PHYLUM ...]]
                          [-2 FASTQ] [--method {flye,unicycler}]
                          [--longread_type {nano-raw,nano-hq,pacbio-raw,pacbio-hifi}]
                          [--genome_size SIZE] [-w DIR] [--tmpdir DIR]
                          [--assembler_args ARG] [-ml BP]
                          [-a [ACCESSION ...]] [-u [URL ...]] [-mc BP]
                          [--sourdb DB ...] [--mincovpct PCT] [-c INT]
                          [-m GB] [-q] [-v]

.. code-block:: bash

    AAFTF pipeline_hybrid \
        -1 reads/STRAINX_R1.fq.gz -2 reads/STRAINX_R2.fq.gz \
        -lr reads/STRAINX_ont.fq.gz \
        -o STRAINX -c 24 -m 96 --phylum Ascomycota

This produces the trim and filter FASTQs as in ``pipeline_short``, then ``STRAINX.flye.fasta``,
``STRAINX.vecscreen.fasta``, ``STRAINX.taxpurge.fasta``, ``STRAINX.covpurge.fasta``,
``STRAINX.rmdup.fasta``, ``STRAINX.racon.fasta``, ``STRAINX.polypolish.fasta``
(``STRAINX.pypolca.fasta`` with single-end reads) and ``STRAINX.final.fasta``. With ``--method unicycler`` the assembly is ``STRAINX.unicycler.fasta``
and no polish files are written.
