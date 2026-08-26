====
trim
====

**Aliases:** ``trim_reads``, ``read_trim``

Adapter and quality trims raw Illumina FASTQ reads (paired- or single-end). This is normally the
first step of the pipeline.

Algorithm
=========

Three interchangeable trimming engines are supported via ``--method``:

* **bbduk** (default)
    Runs BBTools' ``bbduk.sh`` with ``ref=adapters ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo``.
    Because this BBDuk build's paired-mode ``PairStreamer`` (``in1=``/``in2=``) has a bug that
    silently truncates the read stream on large/variable-length paired FASTQ after a few hundred
    reads, AAFTF works around it by first interleaving R1/R2 with ``shuffle.sh`` into one file,
    running BBDuk's (unaffected) single-end reader on the interleaved stream, then
    de-interleaving the trimmed output back to ``_1P``/``_2P`` with ``reformat.sh``.
* **trimmomatic**
    Runs Trimmomatic in ``PE``/``SE`` mode with ``ILLUMINACLIP``, ``LEADING``, ``TRAILING``, and
    ``SLIDINGWINDOW`` steps. AAFTF auto-locates the ``trimmomatic.jar``/adapter files from a
    homebrew shell wrapper or a bioconda Python launcher; pass ``--trimmomatic /path/to.jar`` to
    override.
* **fastp**
    Runs ``fastp`` with ``--low_complexity_filter``, optional 5'/3'/sliding-window quality
    trimming (``--cutfront``/``--cuttail``/``--cutright``), optional PCR-duplicate removal
    (``--dedup``), and optional read merging (``--merge``, writes a ``_MG.fastq.gz`` file of
    merged overlapping pairs). Produces an HTML/JSON fastp QC report alongside the trimmed reads.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Meaning
   * - ``-ml/--minlen``
     - 75
     - Minimum read length retained after trimming
   * - ``-aq/--avgqual``
     - 10
     - Minimum average base quality (bbduk ``maq=``, fastp ``--average_qual``)
   * - ``--method``
     - bbduk
     - bbduk / trimmomatic / fastp
   * - ``-m/--memory``
     - auto (60% of detected system RAM)
     - Max heap for bbduk (``-Xmx``)
   * - Trimmomatic ``LEADING``/``TRAILING``
     - 3 / 3
     - Per-base quality trim from each read end
   * - Trimmomatic ``SLIDINGWINDOW``
     - 4:15
     - Window size:quality for sliding-window trim
   * - Trimmomatic ``ILLUMINACLIP``
     - ``TruSeq3-PE.fa:2:30:10``
     - seed mismatches:palindrome clip threshold:simple clip threshold

Invocation
==========

.. code-block:: text

    AAFTF trim -l LEFT [-r RIGHT] [-o BASENAME] [-c CPUS] [-ml MINLEN] [-aq AVGQUAL]
               [--method {bbduk,trimmomatic,fastp}] [-m MEMORY]
               [--dedup] [--cutfront] [--cuttail] [--cutright] [--merge]
               [--trimmomatic JAR] [--trimmomatic_adaptors FILE] [--trimmomatic_clip STR]
               [--trimmomatic_leadingwindow N] [--trimmomatic_trailingwindow N]
               [--trimmomatic_slidingwindow W:Q] [--trimmomatic_quality {phred33,phred64}]
               [-v] [--pipe]

``-l/--left`` is required; ``-r/--right`` is optional (omit for single-end reads). If
``-o/--out`` (``basename``) is not given, it is derived from ``--left``'s filename (text before
the first ``_`` or ``.``).

**Output:** paired mode writes ``{basename}_1P.fastq.gz`` / ``{basename}_2P.fastq.gz``;
single-end mode writes ``{basename}_1U.fastq.gz``.

Example
=======

.. code-block:: bash

    AAFTF trim --method bbduk --memory 64 -c 16 \
        --left reads/STRAINX_R1.fq.gz --right reads/STRAINX_R2.fq.gz \
        -o reads_trimmed/STRAINX

    # fastp with deduplication and read merging
    AAFTF trim --method fastp -c 16 --dedup --merge \
        --left reads/STRAINX_R1.fq.gz --right reads/STRAINX_R2.fq.gz \
        -o reads_trimmed/STRAINX

Next step: :doc:`filter` (or :doc:`mito` first, for paired data, to seed a mitochondrial
reference used to keep MT reads out of the nuclear filter step).
