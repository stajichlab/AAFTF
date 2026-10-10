======
polish
======

Error-corrects (polishes) the assembly by mapping reads back to it and calling/correcting base
errors, small indels, and local misassemblies.

``polish`` is **optional** and not part of :doc:`pipeline`. It is recommended for assemblies built
from (or supplemented with) long reads; polishing a short-read-only assembly with the same short
reads it was assembled from rarely improves it and can introduce errors.

Algorithm
=========

Four polishing engines are supported via ``--method``, each isolated into its own function
in ``aaftf/polish.py``:

* **polypolish** (default) -- Single-pass short-read polishing via
  `Polypolish <https://github.com/rrwick/Polypolish>`_. AAFTF indexes the assembly with
  ``bwa index``, aligns each read file *separately* with ``bwa mem -a`` (reporting all
  alignments, as Polypolish requires), filters the resulting SAM pairs by insert size with
  ``polypolish filter``, then runs ``polypolish polish``. Requires paired reads
  (``-1/--read1`` and ``-2/--read2``); otherwise an error is raised.
* **pypolca** -- Single-pass polishing via `pypolca <https://github.com/replikation/pypolca>`_,
  a Python reimplementation of MaSuRCA's POLCA algorithm (``bwa`` + ``samtools`` + ``freebayes``,
  run internally by pypolca itself) that works with any modern samtools. AAFTF invokes
  ``pypolca run`` once against the input assembly with both read files and copies out the
  corrected FASTA plus its ``.vcf``/``.report`` outputs (saved as ``{output}.vcf`` /
  ``{output}.pypolca_report.txt``). Requires ``-1/--read1`` (and typically ``-2/--read2``).
* **nextpolish2** -- Repeat-aware polishing of HiFi assemblies via
  `NextPolish2 <https://github.com/Nextomics/NextPolish2>`_. AAFTF builds a short-read k-mer
  database with ``yak count``, maps the HiFi long reads to the assembly with
  ``minimap2 -ax map-hifi`` (piped into ``samtools sort``), then runs ``nextPolish2`` using the
  HiFi alignments plus the k-mer database. Requires both short reads (``-1/--read1`` /
  ``-2/--read2``, used to build the k-mer database) and HiFi long reads (``-lr/--longreads``).
* **racon** -- Long-read polishing via `Racon <https://github.com/isovic/racon>`_. AAFTF maps
  the long reads to the assembly with ``minimap2 -x map-ont`` to produce PAF overlaps, then runs
  ``racon`` using those overlaps to correct the assembly. Requires ``-lr/--longreads``; not
  compatible with short-read-only input.

Required executables per method: ``polypolish`` -- ``bwa``, ``samtools``, ``polypolish``;
``pypolca`` -- ``bwa``, ``samtools``, ``freebayes``, ``pypolca``; ``nextpolish2`` -- ``minimap2``,
``samtools``, ``yak``, ``nextPolish2``; ``racon`` -- ``minimap2``, ``racon``. Missing executables
are detected up front and reported together before any work starts (exit code 2). If the
polisher fails, a ``RuntimeError`` points at ``<workdir>/<method>.log``.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``--method``
     - polypolish
     - polypolish / pypolca / nextpolish2 / racon
   * - ``-m/--memory``
     - 16 (GB)
     - pypolca only: divided by ``-c/--cpus`` for pypolca's per-thread ``-m`` (at least 1G);
       ignored by the other methods

Invocation
==========

.. code-block:: text

    AAFTF polish -i FASTA [-o FASTA] [--method {polypolish,pypolca,nextpolish2,racon}]
                 [-1 FASTQ] [-2 FASTQ] [-lr FASTQ]
                 [-c INT] [-m GB]
                 [-w DIR] [-q] [-v]

``-i/--infile`` is required. ``polypolish`` needs ``-1/--read1`` and ``-2/--read2``; ``pypolca``
needs ``-1/--read1`` (and typically ``-2/--read2``); ``nextpolish2`` needs both short reads and ``-lr/--longreads`` (HiFi);
``racon`` needs ``-lr/--longreads``.

Example
=======

.. code-block:: bash

    # Polypolish (default)
    AAFTF polish -c 24 --memory 96 \
        -i genomes/STRAINX.rmdup.fasta -o genomes/STRAINX.polish.fasta \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz

    # pypolca instead of Polypolish
    AAFTF polish --method pypolca \
        -c 24 --memory 96 \
        -i genomes/STRAINX.rmdup.fasta -o genomes/STRAINX.polish.fasta \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz

    # NextPolish2 (HiFi + short reads for the k-mer database)
    AAFTF polish --method nextpolish2 \
        -c 24 --memory 96 \
        -i genomes/STRAINX.rmdup.fasta -o genomes/STRAINX.polish.fasta \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz \
        --longreads hifi_reads.fastq.gz

    # Racon (ONT/PacBio long reads only)
    AAFTF polish --method racon \
        -c 24 --memory 96 \
        -i genomes/STRAINX.rmdup.fasta -o genomes/STRAINX.polish.fasta \
        --longreads longreads.fastq.gz

Next step: :doc:`sort`.
