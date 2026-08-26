=====
depth
=====

**Aliases:** ``coverage``, ``cov``

Maps reads back to the final assembly and reports per-contig and whole-assembly read-depth
statistics -- both a general sequencing-success/QC summary, and a way to flag contigs whose
coverage is anomalously high (a common signature of a residual contaminant or an
under-represented organelle genome that assembled as a separate, over-covered contig).

Algorithm
=========

1. Counts input reads per FASTQ file (for the report's read-input summary).
2. Maps Illumina reads with ``minimap2 -ax sr`` (default; ``--aligner bwa`` for ``bwa mem``
   instead) and/or long reads with ``minimap2 -ax {map-ont,map-pb,map-hifi}`` (selected via
   ``--longread_type``), sorting each to an indexed BAM. When both read types are supplied, the
   two BAMs are merged with ``samtools merge`` before depth calculation.
3. Runs ``samtools flagstat`` per read type for mapping-rate statistics.
4. Runs ``mosdepth`` on the combined BAM in **quantized mode** (``--quantize``, default
   ``0:1:4:100:200:`` -> bins labeled NO_COVERAGE / LOW_COVERAGE / CALLABLE / HIGH_COVERAGE /
   VERY_HIGH_COVERAGE), producing both a per-contig mean-depth summary and a global coverage
   distribution.
5. Parses ``mosdepth.summary.txt`` for per-contig mean depth and ``mosdepth.global.dist.txt``
   for percent of bases covered at >= 1x.
6. Computes assembly mean depth and its **population** standard deviation across contigs
   (population, not sample, SD is intentional: the contigs *are* the entire assembly, not a
   sample drawn from it -- using sample SD would systematically inflate the outlier threshold).
   Contigs are flagged:

   * **OUTLIER** -- mean depth > assembly mean + **3 SD** (likely contaminant or organellar
     sequence).
   * **ELEVATED** -- mean depth between **2 SD** and 3 SD above the assembly mean (worth manual
     inspection, not automatically flagged as contamination).

   Contigs shorter than ``--min_contig_len`` (default 500 bp) are excluded from the outlier
   statistics (short contigs have noisier depth estimates and would distort the mean/SD).
7. Unless ``--no-plot`` (and matplotlib is available), also renders three plots:
   ``<prefix>.depth_heatmap.<format>`` (per-scaffold coverage-class tiles),
   ``<prefix>.depth_barplot.<format>`` (stacked bar of coverage-class proportions), and
   ``<prefix>.depth_histogram.<format>`` (histogram + boxplot of per-scaffold mean depths),
   in ``--plot-format`` (pdf/svg/png; default pdf).

Report contents
================

The text report (default ``coverage_stats.txt``) has three sections:

1. **Read input summary** -- per-file read counts and ``samtools flagstat`` alignment rates.
2. **Whole-assembly coverage** -- two mean-depth estimates (mosdepth's length-weighted global
   mean, and the unweighted arithmetic mean of per-contig means), plus percent of bases covered
   at >= 1x.
3. **Per-contig depth table**, sorted by depth descending, with OUTLIER/ELEVATED flags as above.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - Elevated threshold
     - mean + 2 SD (hard-coded)
     - Flag contigs worth inspecting
   * - Outlier threshold
     - mean + 3 SD (hard-coded)
     - Flag likely contaminant/organellar contigs
   * - ``--min_contig_len``
     - 500 (bp)
     - Minimum contig length included in outlier statistics
   * - ``--aligner``
     - minimap2
     - minimap2 (default) or bwa, for Illumina reads
   * - ``--illumina_preset``
     - sr
     - minimap2 preset for Illumina reads
   * - ``--longread_type``
     - map-ont
     - minimap2 preset for long reads: map-ont / map-pb / map-hifi
   * - ``--quantize``
     - ``0:1:4:100:200:``
     - mosdepth quantize bin boundaries (colon-separated, trailing colon required)

Invocation
==========

.. code-block:: text

    AAFTF depth -i INPUT [-o OUT] [-l LEFT] [-r RIGHT] [-lr LONGREADS]
               [--aligner {minimap2,bwa}] [--illumina_preset {sr,short}]
               [--longread_type {map-ont,map-pb,map-hifi}]
               [--min_contig_len N] [--quantize STR] [--quantize-labels LABELS]
               [--plot-format {pdf,svg,png}] [--no-plot]
               [-c CPUS] [-w WORKDIR] [-v] [--pipe]

``-i/--input`` is required; provide Illumina reads (``-l``/``-r``), long reads (``-lr``), or both.

Example
=======

.. code-block:: bash

    AAFTF depth -i genome.final.fasta \
        --left reads_1P.fastq.gz --right reads_2P.fastq.gz \
        -c 16 -o coverage_report.txt

    # Illumina + long reads together, PDF plots
    AAFTF depth -i genome.final.fasta \
        --left reads_1P.fastq.gz --right reads_2P.fastq.gz \
        --longreads nanopore.fastq.gz --longread_type map-ont \
        -c 16 -o coverage_report.txt

.. note::
   ``depth`` is not part of ``AAFTF pipeline`` -- run it manually against the final
   (:doc:`sort`-produced) assembly as a QC follow-up, alongside :doc:`assess`.
