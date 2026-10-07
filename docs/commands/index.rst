=================
Command Reference
=================

.. toctree::
   :maxdepth: 1

   dependency
   database
   trim
   mito
   filter
   assemble
   vecscreen
   taxpurge
   covpurge
   fcs_screen
   fcs_gx_purge
   rmdup
   polish
   sort
   assess
   depth
   pipeline
   fix_tbl

.. list-table::
   :header-rows: 1
   :widths: 20 16 46

   * - Canonical name
     - Key tools
     - Purpose
   * - :doc:`dependency`
     - (none -- checks ``$PATH``)
     - Setup: check that AAFTF's external tools and Python packages are installed
   * - :doc:`database`
     - (urllib only)
     - Setup: list or download reference databases into ``$AAFTF_DB``
   * - :doc:`trim`
     - bbduk.sh, trimmomatic, fastp
     - Adapter + quality trimming
   * - :doc:`mito`
     - NOVOPlasty, minimap2
     - (Optional) De novo mitochondrial genome assembly
   * - :doc:`filter`
     - bbduk.sh, bowtie2, bwa, minimap2, samtools
     - Remove contaminant/PhiX reads
   * - :doc:`assemble`
     - spades.py, megahit, unicycler, flye
     - Assemble cleaned reads
   * - :doc:`vecscreen`
     - blastn, makeblastdb
     - Vector/Euk/Prok/Mito contamination screen
   * - :doc:`taxpurge`
     - sourmash
     - Remove contigs classified outside the kept phyla (run by the pipelines only with ``-p``)
   * - :doc:`covpurge`
     - bwa, samtools
     - Remove contigs with low Illumina read coverage
   * - :doc:`fcs_screen`
     - run_fcsadaptor.sh (singularity/docker)
     - (Optional) NCBI FCS-adaptor vector trimming
   * - :doc:`fcs_gx_purge`
     - run_gx.py (NCBI FCS-GX)
     - (Optional) NCBI FCS-GX genomic cross-contamination purge
   * - :doc:`rmdup`
     - minimap2
     - Remove redundant/duplicate contigs
   * - :doc:`polish`
     - polypolish, pypolca, nextPolish2, racon, bwa, samtools, minimap2, yak, freebayes
     - (Optional) Short/long-read assembly polishing
   * - :doc:`sort`
     - (BioPython only)
     - Sort contigs by length, rename headers
   * - :doc:`assess`
     - (BioPython only)
     - Assembly completeness / summary statistics
   * - :doc:`depth`
     - minimap2, bwa, samtools, mosdepth
     - (Optional) Per-contig read-depth report + outlier flags
   * - :doc:`pipeline`
     - those of the steps it runs
     - ``pipeline_short``: trim -> filter -> assemble -> vecscreen -> [taxpurge] -> covpurge -> rmdup -> sort -> assess
       in one command; ``pipeline_long`` (long reads) and ``pipeline_hybrid`` (Illumina + long reads)
   * - :doc:`fix_tbl`
     - (none)
     - Annotation: adjust NCBI ``.tbl`` feature coordinates after FCS trimming

``AAFTF -h`` lists the subcommands in four groups: *Setup* (``dependency``, ``database``),
*Steps* (``trim`` through ``depth``, in the order above; steps marked "(Optional)" are not run by
``pipeline_short``), *Pipeline* (``pipeline_short``, ``pipeline_long`` and ``pipeline_hybrid``) and
*Annotation* (``fix_tbl``).

Common options and exit codes
=============================

Every subcommand accepts ``-v/--verbose`` (verbose logging; also usually retains temp working
directories) and ``-q/--quiet`` (only warnings and errors, which also hides the "your next
command might be" hint printed at the end of a step).

Each step with a working directory also writes its log to ``<workdir>/<subcommand>.log`` (everything
at INFO and above, even with ``-q``, plus full tracebacks on errors). The log is kept whenever the
working directory is: when you pass your own ``-w/--workdir``, with ``-v``, or when the step fails;
an auto-created working directory and its log are removed after a successful run. Subcommands
without a working directory (e.g. ``trim``, ``sort``, ``assess``) write ``./<subcommand>.log`` only
with ``-v``.

``AAFTF`` exits with status 0 on success, 1 on an error, 2 when a required input file, database or
tool is missing, and 130 when interrupted with Ctrl-C.
