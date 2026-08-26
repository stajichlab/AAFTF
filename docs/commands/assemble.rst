========
assemble
========

**Aliases:** ``asm``, ``spades``

Runs a de novo genome assembler on cleaned (trimmed + filtered) reads.

Algorithm
=========

The assembler is selected with ``--method``:

* **spades** (default) -- ``spades.py`` with ``--isolate`` and/or ``--careful`` mode (both
  default on; disable with ``--no-isolate``/``--no-careful``), ``--cov-cutoff auto`` (unless
  ``--meta`` is passed via ``--assembler_args``), and optional merged/long reads. If the
  ``--workdir`` already exists from a prior run, AAFTF instead resumes with
  ``spades.py --restart-from last`` rather than restarting from scratch. Final assembly is copied
  from SPAdes' ``scaffolds.fasta``.
* **dipspades** -- SPAdes' diploid-aware assembler (only available in SPAdes <= 3.11.1; not
  packaged in later SPAdes releases). Supports resuming via ``dipspades.py --continue``. Final
  assembly is copied from ``consensus_contigs.fasta``, with paired/unpaired haplotig FASTAs also
  copied out if ``--haplocontigs`` matched files under the ``dipspades/`` subdirectory.
* **megahit** -- fast De Bruijn graph assembler, generally lower accuracy than SPAdes but much
  faster/lower memory; does not support resuming (an existing workdir with the same name errors).
  Final assembly is copied from ``final.contigs.fa``.
* **unicycler** -- wraps SPAdes with additional scaffolding logic; supports combining
  short reads with ``--longreads`` (hybrid assembly). Final assembly is copied from
  ``assembly.fasta``.
* **masurca** / **nextdenovo** -- accepted as ``--method`` values but not yet implemented; AAFTF
  will print a message and exit without assembling.

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
     - spades / dipspades / megahit / unicycler
   * - ``-m/--memory``
     - 32 (GB)
     - Passed to SPAdes ``--mem`` / megahit ``--memory``
   * - ``--careful``
     - on
     - SPAdes ``--careful`` (more accurate, slower; mismatch/indel correction)
   * - ``--isolate``
     - off unless explicitly set*
     - SPAdes ``--isolate`` mode (recommended for high-coverage, low-diversity isolate data)
   * - SPAdes coverage cutoff
     - ``auto``
     - ``--cov-cutoff auto`` unless running in ``--meta`` mode

\* ``--isolate``/``--no-isolate`` are both ``store_true``/``store_false`` onto the same
``isolate`` destination with no shared default set at the parser level; check
``AAFTF assemble --help`` for the resolved default in your installed version, or pass one flag
explicitly to be certain.

Invocation
==========

.. code-block:: text

    AAFTF assemble --method {spades,dipspades,megahit,unicycler} -o OUT
                   [-w WORKDIR] [-c CPUS] [-m MEMORY]
                   [-l LEFT] [-r RIGHT] [-lr LONGREADS] [--single/--merged FILE]
                   [--careful/--no-careful] [--isolate/--no-isolate]
                   [--tmpdir DIR] [--assembler_args ARG ...] [--haplocontigs FILE]
                   [-v] [--pipe]

``-o/--out`` is required (output assembly FASTA path). ``--assembler_args`` may be repeated to
pass through additional raw SPAdes/megahit arguments.

Example
=======

.. code-block:: bash

    AAFTF assemble -c 24 --memory 96 \
        --left reads_filtered_1.fastq.gz --right reads_filtered_2.fastq.gz \
        -o genomes/STRAINX.spades.fasta -w working_AAFTF/spades_STRAINX

    # Hybrid short+long read assembly with Unicycler
    AAFTF assemble --method unicycler -c 24 \
        --left reads_filtered_1.fastq.gz --right reads_filtered_2.fastq.gz \
        --longreads ont_reads.fastq.gz \
        -o genomes/STRAINX.unicycler.fasta

Next step: :doc:`vecscreen`.
