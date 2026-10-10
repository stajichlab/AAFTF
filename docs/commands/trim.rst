====
trim
====

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
    ``SLIDINGWINDOW`` steps. AAFTF auto-locates ``trimmomatic.jar`` from the ``trimmomatic``
    launcher on ``$PATH`` (homebrew shell wrapper or bioconda Python launcher), and looks for the
    adaptor file (``--trimmomatic_adaptors``) next to the jar and under ``<prefix>/share/trimmomatic``
    in the directories above it; a missing jar or adaptor file raises an error (exit code 2).
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
     - 8 (GB)
     - Max Java heap for bbduk (``-Xmx``)
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

    AAFTF trim -1 FASTQ [-2 FASTQ] [-o PREFIX] [-c INT] [-ml BP] [-aq INT]
               [--method {bbduk,trimmomatic,fastp}] [-m GB]
               [--dedup] [--cutfront] [--cuttail] [--cutright] [--merge]
               [--trimmomatic_adaptors FILE] [--trimmomatic_clip STR]
               [--trimmomatic_leadingwindow N] [--trimmomatic_trailingwindow N]
               [--trimmomatic_slidingwindow W:Q] [--trimmomatic_quality {phred33,phred64}]
               [-q] [-v]

``-1/--read1`` is required; ``-2/--read2`` is optional (omit for single-end reads). If
``-o/--out`` (``basename``) is not given, it is derived from ``--read1``'s filename (text before
the first ``_``, or the first ``.`` if it has none). Supplying no reads raises an error.

**Output:** paired mode writes ``{basename}_1P.fastq.gz`` / ``{basename}_2P.fastq.gz``;
single-end mode writes ``{basename}_1U.fastq.gz``.

Example
=======

.. code-block:: bash

    AAFTF trim --method bbduk --memory 64 -c 16 \
        --read1 reads/STRAINX_R1.fq.gz --read2 reads/STRAINX_R2.fq.gz \
        -o reads_trimmed/STRAINX

    # fastp with deduplication and read merging
    AAFTF trim --method fastp -c 16 --dedup --merge \
        --read1 reads/STRAINX_R1.fq.gz --read2 reads/STRAINX_R2.fq.gz \
        -o reads_trimmed/STRAINX

Next step: :doc:`filter` (or :doc:`mito` first, for paired data, to seed a mitochondrial
reference used to keep MT reads out of the nuclear filter step).
