=================
Command Reference
=================

.. toctree::
   :maxdepth: 1

   download
   trim
   mito
   filter
   assemble
   vecscreen
   fcs_screen
   fcs_gx_purge
   sourpurge
   rmdup
   polish
   sort
   assess
   depth
   fix_tbl
   pipeline

.. list-table::
   :header-rows: 1
   :widths: 16 22 16 46

   * - Canonical name
     - Aliases
     - Key tools
     - Purpose
   * - :doc:`download`
     - ``configure``, ``install``, ``download_db``, ``setup``
     - (urllib only)
     - Fetch/cache reference databases into ``$AAFTF_DB``
   * - :doc:`trim`
     - ``trim_reads``, ``read_trim``
     - bbduk.sh, trimmomatic, fastp
     - Adapter + quality trimming
   * - :doc:`mito`
     - ``mito_asm``, ``mitochondria``
     - NOVOPlasty, minimap2
     - De novo mitochondrial genome assembly
   * - :doc:`filter`
     - ``filter_reads``, ``read_filter``
     - bbduk.sh, bowtie2, bwa, minimap2, samtools
     - Remove contaminant/PhiX reads
   * - :doc:`assemble`
     - ``asm``, ``spades``
     - spades.py, megahit, unicycler
     - Assemble cleaned reads
   * - :doc:`vecscreen`
     - ``vectorscreen``, ``vector_blast``
     - blastn, makeblastdb
     - Vector/Euk/Prok/Mito contamination screen
   * - :doc:`fcs_screen`
     - ``ncbi_fcs``, ``ncbi_fcs-screen``
     - run_fcsadaptor.sh (singularity/docker)
     - NCBI FCS-adaptor vector trimming
   * - :doc:`fcs_gx_purge`
     - ``ncbi_fcs-gx``, ``ncbi_fcs_gx``, ``gx``
     - run_gx.py (NCBI FCS-GX)
     - NCBI FCS-GX genomic cross-contamination purge
   * - :doc:`sourpurge`
     - ``purge``
     - sourmash, bwa, samtools
     - Taxonomy + low-coverage contig purge
   * - :doc:`rmdup`
     - ``dedup``
     - minimap2
     - Remove redundant/duplicate contigs
   * - :doc:`polish`
     - ``pilon``, ``polca``
     - pilon, nextPolish, polca, bwa, samtools
     - Short/long-read assembly polishing
   * - :doc:`sort`
     - --
     - (BioPython only)
     - Sort contigs by length, rename headers
   * - :doc:`assess`
     - ``stats``
     - (BioPython only)
     - Assembly completeness / summary statistics
   * - :doc:`depth`
     - ``coverage``, ``cov``
     - minimap2, bwa, samtools, mosdepth
     - Per-contig read-depth report + outlier flags
   * - :doc:`fix_tbl`
     - ``fix``
     - (none)
     - Adjust NCBI ``.tbl`` feature coordinates after FCS trimming
   * - :doc:`pipeline`
     - --
     - all of the above (except fcs_screen/fcs_gx_purge/depth)
     - Run the full trim -> ... -> assess pipeline in one command

Every subcommand accepts ``-v/--debug`` (verbose logging; also usually retains temp working
directories) and ``--pipe`` (suppress the "your next command might be" hint printed at the end --
set automatically when a step is invoked from inside ``AAFTF pipeline``).
