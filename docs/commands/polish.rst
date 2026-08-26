======
polish
======

**Aliases:** ``pilon``, ``polca``

Error-corrects (polishes) the assembly by mapping reads back to it and calling/correcting base
errors, small indels, and local misassemblies.

Algorithm
=========

Four interchangeable polishing engines are supported via ``--method``:

* **pilon** (default) -- Iterative, short-read polishing. Each iteration: build a fresh
  ``bwa mem`` + ``samtools sort``/``markdup`` BAM of reads against the current assembly (keeping
  only properly-paired reads, SAM flag ``0x2``), then run Pilon with ``--frags BAMFILE
  --changes``. The number of changes Pilon reports is used as an automatic early-stopping signal:
  if an iteration makes **zero** changes, polishing stops early even if
  ``-it/--iterations`` has not been reached. ``--diploid`` (or ``--ploidy 2``) switches Pilon to
  diploid heterozygous-SNP handling.
* **nextpolish** -- Two-task NextPolish correction per iteration (task 1 then task 2, each
  against a freshly-remapped BAM), also honoring ``--ploidy``/``--diploid``.
* **polca** / **masurca** -- Single-pass polishing via MaSuRCA's ``polca.sh``, which internally
  handles its own read mapping and iterative correction; AAFTF just invokes it once against
  the input assembly with both read files and copies out the corrected FASTA plus its
  ``.vcf``/``.report`` outputs (saved as ``{output}.vcf`` / ``{output}.polca_report.txt``).
* **racon** -- Long-read polishing mode (requires ``-lr/--longreads``; not compatible with
  short-read-only input).

Every short-read method requires ``bwa`` and ``samtools`` on ``$PATH`` in addition to the chosen
polisher itself; missing executables are detected up front and reported together before any work
starts.

.. warning::
   **polca.sh + modern samtools.** MaSuRCA's bundled ``polca.sh`` uses ``samtools sort -f``,
   which was removed in samtools 1.21+. If system samtools is >= 1.21, AAFTF prints a warning
   (and, on failure, an error) suggesting either ``--polca_samtools /path/to/older/samtools``
   (points polca specifically at a compatible binary via the ``SAMTOOLS`` env var and a PATH
   prefix) or switching to ``--method pilon``. The project also ships a pre-patched
   ``patches/polca.sh`` (already applied automatically in the Docker/Singularity images) that
   fixes this at the source.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``--method``
     - pilon
     - pilon / polca / nextpolish / racon
   * - ``-it/--iterations``
     - 5
     - Max Pilon/NextPolish rounds (pilon may stop earlier once 0 changes are reported)
   * - ``-m/--memory``
     - 16 (GB)
     - Total memory budget; divided by ``-c/--cpus`` for per-thread BAM-sort memory
   * - ``--ploidy``
     - 1
     - NextPolish ploidy (``--diploid``/``--ploidy 2`` forces diploid mode for Pilon too)

Invocation
==========

.. code-block:: text

    AAFTF polish -i INFILE [-o OUTFILE] --method {pilon,polca,nextpolish,racon}
                [-l LEFT] [-r RIGHT] [-lr LONGREADS]
                [-c CPUS] [-m MEMORY] [-it ITERATIONS]
                [--diploid] [--ploidy N] [--polca PATH] [--polca_samtools PATH]
                [-w WORKDIR] [-v] [--pipe]

``-i/--infile`` is required. Short-read methods (pilon/polca/nextpolish) need ``-l/--left``
(and typically ``-r/--right``); ``racon`` needs ``-lr/--longreads``.

Example
=======

.. code-block:: bash

    AAFTF polish -c 24 --memory 96 -it 5 \
        -i genomes/STRAINX.rmdup.fasta -o genomes/STRAINX.polish.fasta \
        --left reads_filtered_1.fastq.gz --right reads_filtered_2.fastq.gz

    # POLCA instead of Pilon, pointing at an older samtools for compatibility
    AAFTF polish --method polca --polca_samtools /opt/samtools-1.17/bin/samtools \
        -c 24 --memory 96 \
        -i genomes/STRAINX.rmdup.fasta -o genomes/STRAINX.polish.fasta \
        --left reads_filtered_1.fastq.gz --right reads_filtered_2.fastq.gz

Next step: :doc:`sort`.
