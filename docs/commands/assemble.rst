========
assemble
========

Runs a de novo genome assembler on cleaned (trimmed + filtered) short reads, or on long reads
(Flye).

Algorithm
=========

The assembler is selected with ``--method``, which is required (there is no default assembler;
``pipeline_short`` uses spades, ``pipeline_long`` and ``pipeline_hybrid`` use flye):

* **spades** -- ``spades.py --mem <GB>`` in ``--isolate`` mode (default; disable with
  ``--no-isolate``). SPAdes does not allow ``--careful`` together with ``--isolate``, so
  ``--careful`` is only added with ``--no-isolate`` (and dropped with ``--no-careful``).
  ``--cov-cutoff auto`` is added unless ``--meta`` is passed via ``--assembler_args``; merged/single
  reads from ``--merged`` are included. If the
  ``--workdir`` already exists from a prior run, AAFTF instead resumes with
  ``spades.py --restart-from last`` rather than restarting from scratch. Final assembly is copied
  from SPAdes' ``scaffolds.fasta``.
* **megahit** -- fast De Bruijn graph assembler, generally lower accuracy than SPAdes but much
  faster/lower memory; ``-m`` GB is converted to bytes for ``megahit --memory``. Does not support
  resuming: if the output folder already exists, AAFTF stops with an error before running MEGAHIT
  (remove it or pass another ``-w``). MEGAHIT has no merged-read option, so ``--merged`` reads
  are passed with ``-r`` as single-end reads (alongside ``-1``/``-2``, or added to the
  single-end ``-r`` list).
  Final assembly is copied from ``final.contigs.fa``.
* **unicycler** -- wraps SPAdes with additional scaffolding logic; supports combining
  short reads with ``--longreads`` (hybrid assembly). With paired reads, ``--merged`` reads are
  passed as ``--unpaired`` alongside the pairs (ignored for single-end input). Final assembly is
  copied from ``assembly.fasta``.
* **flye** -- long-read assembler for ONT/PacBio reads (``-lr/--longreads``, required; ``-1`` is
  not used). ``--longread_type`` selects Flye's ``--nano-raw`` / ``--nano-hq`` (default) /
  ``--pacbio-raw`` / ``--pacbio-hifi``; ``--genome_size`` (e.g. ``40m``) is passed as
  ``--genome-size`` when given. If the ``--workdir`` holds a previous Flye run (``params.json``),
  Flye is restarted with ``--resume``. Final assembly is copied from ``assembly.fasta``; without
  ``-o`` it is named ``<long-read prefix>.flye.fasta``. The next-step hint suggests
  ``AAFTF polish --method racon``.

If the assembler's expected output file is missing after the run, a ``RuntimeError`` is raised
(exit code 1) pointing at the assembler log in the work directory.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``--method``
     - spades
     - spades / megahit / unicycler / flye
   * - ``-m/--memory``
     - 32 (GB)
     - Integer GB; SPAdes ``--mem`` / MEGAHIT ``--memory`` (converted to bytes)
   * - ``--no-isolate``
     - isolate on
     - SPAdes ``--isolate`` mode (recommended for high-coverage isolate data)
   * - ``--no-careful``
     - careful on
     - SPAdes ``--careful``; only applied when ``--no-isolate`` is given
   * - SPAdes coverage cutoff
     - ``auto``
     - ``--cov-cutoff auto`` unless running in ``--meta`` mode
   * - ``--longread_type``
     - nano-hq
     - Flye read type: nano-raw / nano-hq / pacbio-raw / pacbio-hifi
   * - ``--genome_size``
     - estimated by Flye
     - Flye ``--genome-size``, e.g. ``40m``

Invocation
==========

.. code-block:: text

    AAFTF assemble -o FASTA --method {spades,megahit,unicycler,flye} [-1 FASTQ] [-2 FASTQ]
                   [-w DIR] [-c INT] [-m GB] [--merged FASTQ] [-lr FASTQ]
                   [--longread_type {nano-raw,nano-hq,pacbio-raw,pacbio-hifi}]
                   [--genome_size SIZE]
                   [--no-careful] [--no-isolate]
                   [--tmpdir DIR] [--assembler_args ARG]
                   [-q] [-v]

``-o/--out`` (output assembly FASTA) and ``--method`` are required, and ``-1/--read1`` for spades,
megahit and unicycler. ``--assembler_args`` takes
one argument and may be repeated to pass several raw arguments to the assembler. ``-lr/--longreads``
is required for Flye and optional for Unicycler. The step log is written to ``<workdir>/assemble.log``.

Example
=======

.. code-block:: bash

    AAFTF assemble --method spades -c 24 --memory 96 \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz \
        -o genomes/STRAINX.spades.fasta -w working_AAFTF/spades_STRAINX

    # Hybrid short+long read assembly with Unicycler
    AAFTF assemble --method unicycler -c 24 \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz \
        --longreads ont_reads.fastq.gz \
        -o genomes/STRAINX.unicycler.fasta

    # Long-read assembly with Flye (ONT reads, ~40 Mb genome)
    AAFTF assemble --method flye -c 24 \
        --longreads ont_reads.fastq.gz --genome_size 40m \
        -o genomes/STRAINX.flye.fasta -w working_AAFTF/flye_STRAINX

Next step: :doc:`vecscreen` (for long-read assemblies, :doc:`polish` ``--method racon`` after :doc:`rmdup`).
