========================================
Overall Workflow: Reads to Clean Genome
========================================

AAFTF chains a fixed sequence of subcommands to go from raw Illumina (optionally + long-read)
sequencing data to a polished, sorted, and QC'd genome assembly. Each step consumes the previous
step's output and, when not run with ``--pipe``, prints the exact next AAFTF command to run --
so the whole pipeline can be worked through interactively one command at a time, or run end-to-end
with the single ``AAFTF pipeline`` subcommand (see :doc:`commands/pipeline`).

Pipeline diagram
=================

.. code-block:: text

    raw FASTQ (R1/R2, + optional long reads)
        |
        v
    [1] trim ---------------------------- BBDuk / fastp / Trimmomatic adapter+quality trim
        |
        +--> [1b] mito (optional, PE only) -- de novo NOVOPlasty mitochondrial assembly
        |                                      (used as an extra screen_local reference below)
        v
    [2] filter --------------------------- remove PhiX/UniVec/user-specified contaminant reads
        |
        v
    [3] assemble -------------------------- SPAdes (default) / megahit / unicycler / dipSPAdes
        |
        v
    [4] vecscreen -------------------------- BLASTN vector + Euk/Prok/Mito contamination screen
        |
        v
    [5a] sourpurge   (or)  [5b] fcs_gx_purge / fcs_screen
        sourmash LCA taxonomy    NCBI Foreign Contamination Screen (adaptor and/or genomic-cross
        + low-coverage purge     contamination via NCBI's fcs-adaptor / fcs-gx tools)
        |
        v
    [6] rmdup ------------------------------ minimap2 self-alignment removes redundant contigs
        |
        v
    [7] polish ------------------------------ Pilon / POLCA / NextPolish / Racon short/long-read
        |                                      error correction (iterative)
        v
    [8] sort -------------------------------- rank contigs longest->shortest, rename headers
        |
        v
    [9] assess ------------------------------- N50/L50/GC%/telomere/gap statistics report
        |
        v
    [10] depth (optional) -------------------- per-contig read-depth report + outlier/contaminant
                                                 flags, coverage plots

Expected inputs/outputs per step
==================================

.. list-table::
   :header-rows: 1
   :widths: 14 30 30

   * - Step
     - Input
     - Output
   * - trim
     - raw FASTQ
     - ``{base}_1P.fastq.gz``, ``{base}_2P.fastq.gz``
   * - mito (optional)
     - trimmed PE FASTQ
     - ``{base}.mito.fasta``
   * - filter
     - trimmed FASTQ
     - ``{base}_filtered_1.fastq.gz``, ``{base}_filtered_2.fastq.gz``
   * - assemble
     - filtered FASTQ
     - ``{base}.{method}.fasta``
   * - vecscreen
     - assembled FASTA
     - ``{base}.vecscreen.fasta`` (+ ``{base}.mitochondria.fasta``)
   * - sourpurge / fcs_gx_purge
     - vecscreen FASTA
     - ``{base}.sourpurge.fasta``
   * - rmdup
     - sourpurge FASTA
     - ``{base}.rmdup.fasta``
   * - polish
     - rmdup FASTA
     - ``{base}.polish.fasta``
   * - sort
     - polished FASTA
     - ``{base}.final.fasta``
   * - assess
     - sorted FASTA
     - printed stats (optional ``--report`` file)
   * - depth
     - final FASTA + reads
     - ``coverage_stats.txt`` (+ plots)

Running the whole thing end-to-end
=====================================

The :doc:`commands/pipeline` subcommand runs trim -> mito (if paired) -> filter -> assemble ->
vecscreen -> sourpurge -> rmdup -> polish -> sort -> assess automatically, skipping any step whose
output file already exists (so an interrupted run can simply be re-launched):

.. code-block:: bash

    AAFTF pipeline \
        -l reads_R1.fq.gz -r reads_R2.fq.gz \
        -o STRAINX -c 24 -m 96 \
        --phylum Ascomycota \
        --AAFTF_DB "$AAFTF_DB"

``fcs_screen``/``fcs_gx_purge`` and ``depth`` are *not* part of ``pipeline`` -- run them
separately (see :doc:`commands/fcs_screen`, :doc:`commands/fcs_gx_purge`,
:doc:`commands/depth`) if you want NCBI FCS-based contamination screening in place of/alongside
sourmash, or a coverage report of the final assembly.

Running step-by-step
=======================

.. code-block:: bash

    MEM=128 CPU=24 BASE=STRAINX
    READSDIR=reads TRIMREAD=reads_trimmed OUTDIR=genomes
    mkdir -p "$TRIMREAD" "$OUTDIR"

    AAFTF trim --method bbduk --memory $MEM -c $CPU \
        --left $READSDIR/${BASE}_R1.fq.gz --right $READSDIR/${BASE}_R2.fq.gz \
        -o $TRIMREAD/${BASE}

    AAFTF filter -c $CPU --memory $MEM --aligner bbduk \
        -o $TRIMREAD/${BASE} \
        --left $TRIMREAD/${BASE}_1P.fastq.gz --right $TRIMREAD/${BASE}_2P.fastq.gz

    AAFTF assemble -c $CPU --memory $MEM \
        --left $TRIMREAD/${BASE}_filtered_1.fastq.gz --right $TRIMREAD/${BASE}_filtered_2.fastq.gz \
        -o $OUTDIR/${BASE}.spades.fasta -w working_AAFTF/spades_${BASE}

    AAFTF vecscreen -c $CPU -i $OUTDIR/${BASE}.spades.fasta -o $OUTDIR/${BASE}.vecscreen.fasta

    AAFTF sourpurge -c $CPU --phylum Ascomycota \
        -i $OUTDIR/${BASE}.vecscreen.fasta -o $OUTDIR/${BASE}.sourpurge.fasta \
        --left $TRIMREAD/${BASE}_filtered_1.fastq.gz --right $TRIMREAD/${BASE}_filtered_2.fastq.gz

    AAFTF rmdup -c $CPU -i $OUTDIR/${BASE}.sourpurge.fasta -o $OUTDIR/${BASE}.rmdup.fasta

    AAFTF polish -c $CPU --memory $MEM -i $OUTDIR/${BASE}.rmdup.fasta -o $OUTDIR/${BASE}.polish.fasta \
        --left $TRIMREAD/${BASE}_filtered_1.fastq.gz --right $TRIMREAD/${BASE}_filtered_2.fastq.gz

    AAFTF sort -i $OUTDIR/${BASE}.polish.fasta -o $OUTDIR/${BASE}.final.fasta

    AAFTF assess -i $OUTDIR/${BASE}.final.fasta -r $OUTDIR/${BASE}.stats.txt

    AAFTF depth -i $OUTDIR/${BASE}.final.fasta \
        --left $TRIMREAD/${BASE}_filtered_1.fastq.gz --right $TRIMREAD/${BASE}_filtered_2.fastq.gz \
        -c $CPU -o $OUTDIR/${BASE}.coverage_stats.txt

Each individual command page under :doc:`commands/index` documents its own defaults and cutoffs
in detail.
