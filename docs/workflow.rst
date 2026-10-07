========================================
Overall Workflow: Reads to Clean Genome
========================================

AAFTF chains a fixed sequence of subcommands to go from raw sequencing reads to a cleaned, sorted
and QC'd genome assembly. There are three workflows, one per kind of data, each with its own
all-in-one subcommand (see :doc:`commands/pipeline`):

* :ref:`workflow-short` -- Illumina reads only: ``pipeline_short``
* :ref:`workflow-long` -- long reads (ONT/PacBio) only: ``pipeline_long``
* :ref:`workflow-hybrid` -- Illumina + long reads: ``pipeline_hybrid``

Each step consumes the previous step's output and, unless run with ``-q/--quiet``, prints the
exact next AAFTF command to run, so a workflow can be worked through one command at a time or run
end-to-end with its ``pipeline_*`` subcommand. A pipeline skips any step whose output file already
exists, so an interrupted run can simply be re-launched with the same command.

All three workflows clean the contigs before anything else is done with them: :doc:`commands/vecscreen`
(vector and contaminant BLAST screen), :doc:`commands/taxpurge` (sourmash taxonomy screen, run only
with ``-p/--phylum``; leave it out if the organism's taxonomy is unknown), :doc:`commands/covpurge`
(low Illumina read coverage; not in the long-read workflow) and :doc:`commands/rmdup` (redundant
contigs). Polishing, when the workflow has it, comes after cleaning so discarded contigs are never
polished.

In the diagrams below, dashed boxes are optional steps and dashed arrows are the alternative paths
(e.g. without ``--phylum``, taxpurge is skipped).

.. _workflow-short:

Short reads (pipeline_short)
=============================

.. mermaid::

    flowchart TD
        reads[/"Illumina FASTQ<br/>R1/R2 or single-end"/]
        trim["<b>trim</b><br/>BBDuk / fastp / Trimmomatic<br/>adapter + quality trim"]
        mito["<b>mito</b> (optional, PE only)<br/>NOVOPlasty mitochondrial assembly"]:::opt
        filter["<b>filter</b><br/>remove PhiX / UniVec / user<br/>contaminant reads"]
        assemble["<b>assemble</b><br/>spades (default) / megahit / unicycler"]
        vecscreen["<b>vecscreen</b><br/>BLASTN vector + Euk/Prok/Mito screen"]
        taxpurge["<b>taxpurge</b> (with --phylum)<br/>sourmash taxonomy screen"]:::opt
        fcs["<b>fcs_screen / fcs_gx_purge</b><br/>NCBI FCS (run separately)"]:::opt
        covpurge["<b>covpurge</b><br/>drop low-coverage contigs"]
        rmdup["<b>rmdup</b><br/>remove redundant contigs"]
        sort["<b>sort</b><br/>longest first, rename headers"]
        assess["<b>assess</b><br/>N50 / GC / telomere / gap stats"]
        depth["<b>depth</b> (optional)<br/>per-contig read depth"]:::opt

        reads --> trim --> filter --> assemble --> vecscreen --> taxpurge --> covpurge --> rmdup --> sort --> assess --> depth
        trim -.-> mito -. "--screen_local" .-> filter
        vecscreen -. "no --phylum" .-> covpurge
        vecscreen -.-> fcs -.-> covpurge
    classDef opt stroke-dasharray: 5 5

Polishing a short-read-only assembly with the same reads rarely helps, so ``pipeline_short`` does
not run :doc:`commands/polish`.

.. list-table:: Files written (``{base}`` is the ``-o`` prefix)
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
   * - taxpurge (with ``--phylum``)
     - vecscreen FASTA
     - ``{base}.taxpurge.fasta``
   * - covpurge
     - taxpurge (or vecscreen) FASTA
     - ``{base}.covpurge.fasta``
   * - rmdup
     - covpurge FASTA
     - ``{base}.rmdup.fasta``
   * - sort
     - rmdup FASTA
     - ``{base}.final.fasta``
   * - assess
     - sorted FASTA
     - printed stats (optional ``--report`` file)
   * - depth (optional)
     - final FASTA + reads
     - ``coverage_stats.txt`` (+ plots)

One command:

.. code-block:: bash

    AAFTF pipeline_short \
        -1 reads_R1.fq.gz -2 reads_R2.fq.gz \
        -o STRAINX -c 24 -m 96 \
        --phylum Ascomycota

Step by step:

.. code-block:: bash

    MEM=128 CPU=24 BASE=STRAINX
    READSDIR=reads TRIMREAD=reads_trimmed OUTDIR=genomes
    mkdir -p "$TRIMREAD" "$OUTDIR"

    AAFTF trim --method bbduk --memory $MEM -c $CPU \
        --read1 $READSDIR/${BASE}_R1.fq.gz --read2 $READSDIR/${BASE}_R2.fq.gz \
        -o $TRIMREAD/${BASE}

    AAFTF filter -c $CPU --memory $MEM --aligner bbduk \
        -o $TRIMREAD/${BASE} \
        --read1 $TRIMREAD/${BASE}_1P.fastq.gz --read2 $TRIMREAD/${BASE}_2P.fastq.gz

    AAFTF assemble --method spades -c $CPU --memory $MEM \
        --read1 $TRIMREAD/${BASE}_filtered_1.fastq.gz --read2 $TRIMREAD/${BASE}_filtered_2.fastq.gz \
        -o $OUTDIR/${BASE}.spades.fasta -w working_AAFTF/spades_${BASE}

    AAFTF vecscreen -c $CPU -i $OUTDIR/${BASE}.spades.fasta -o $OUTDIR/${BASE}.vecscreen.fasta

    # skip taxpurge if the organism's phylum is unknown (and pass the vecscreen FASTA to covpurge)
    AAFTF taxpurge --phylum Ascomycota \
        -i $OUTDIR/${BASE}.vecscreen.fasta -o $OUTDIR/${BASE}.taxpurge.fasta

    AAFTF covpurge -c $CPU \
        -i $OUTDIR/${BASE}.taxpurge.fasta -o $OUTDIR/${BASE}.covpurge.fasta \
        --read1 $TRIMREAD/${BASE}_filtered_1.fastq.gz --read2 $TRIMREAD/${BASE}_filtered_2.fastq.gz

    AAFTF rmdup -c $CPU -i $OUTDIR/${BASE}.covpurge.fasta -o $OUTDIR/${BASE}.rmdup.fasta

    AAFTF sort -i $OUTDIR/${BASE}.rmdup.fasta -o $OUTDIR/${BASE}.final.fasta

    AAFTF assess -i $OUTDIR/${BASE}.final.fasta -r $OUTDIR/${BASE}.stats.txt

    # optional coverage report of the final assembly
    AAFTF depth -i $OUTDIR/${BASE}.final.fasta \
        --read1 $TRIMREAD/${BASE}_filtered_1.fastq.gz --read2 $TRIMREAD/${BASE}_filtered_2.fastq.gz \
        -c $CPU -o $OUTDIR/${BASE}.coverage_stats.txt

.. _workflow-long:

Long reads (pipeline_long)
===========================

.. mermaid::

    flowchart TD
        reads[/"long-read FASTQ<br/>ONT or PacBio"/]
        assemble["<b>assemble</b><br/>flye (--longread_type, --genome_size)"]
        vecscreen["<b>vecscreen</b><br/>BLASTN vector + Euk/Prok/Mito screen"]
        taxpurge["<b>taxpurge</b> (with --phylum)<br/>sourmash taxonomy screen"]:::opt
        rmdup["<b>rmdup</b><br/>remove redundant contigs"]
        polish["<b>polish</b><br/>racon with the long reads"]
        sort["<b>sort</b><br/>longest first, rename headers"]
        assess["<b>assess</b><br/>N50 / GC / telomere / gap stats"]
        depth["<b>depth</b> (optional)<br/>per-contig read depth"]:::opt

        reads --> assemble --> vecscreen --> taxpurge --> rmdup --> polish --> sort --> assess --> depth
        vecscreen -. "no --phylum" .-> rmdup
        reads -. "long reads" .-> polish
    classDef opt stroke-dasharray: 5 5

There are no reads to trim or filter, and :doc:`commands/covpurge` is not run, since it needs
Illumina reads for contig coverage.

.. list-table:: Files written (``{base}`` is the ``-o`` prefix)
   :header-rows: 1
   :widths: 14 30 30

   * - Step
     - Input
     - Output
   * - assemble
     - long-read FASTQ
     - ``{base}.flye.fasta``
   * - vecscreen
     - flye FASTA
     - ``{base}.vecscreen.fasta`` (+ ``{base}.mitochondria.fasta``)
   * - taxpurge (with ``--phylum``)
     - vecscreen FASTA
     - ``{base}.taxpurge.fasta``
   * - rmdup
     - taxpurge (or vecscreen) FASTA
     - ``{base}.rmdup.fasta``
   * - polish (racon)
     - rmdup FASTA + long reads
     - ``{base}.racon.fasta``
   * - sort
     - racon FASTA
     - ``{base}.final.fasta``
   * - assess
     - sorted FASTA
     - printed stats (optional ``--report`` file)

One command:

.. code-block:: bash

    AAFTF pipeline_long \
        -lr reads/STRAINX_ont.fq.gz \
        -o STRAINX -c 24 --genome_size 40m \
        --phylum Ascomycota

Step by step:

.. code-block:: bash

    CPU=24 BASE=STRAINX LONGREADS=reads/STRAINX_ont.fq.gz OUTDIR=genomes
    mkdir -p "$OUTDIR"

    AAFTF assemble --method flye -c $CPU \
        --longreads $LONGREADS --longread_type nano-hq --genome_size 40m \
        -o $OUTDIR/${BASE}.flye.fasta -w working_AAFTF/flye_${BASE}

    AAFTF vecscreen -c $CPU -i $OUTDIR/${BASE}.flye.fasta -o $OUTDIR/${BASE}.vecscreen.fasta

    # skip taxpurge if the organism's phylum is unknown (and pass the vecscreen FASTA to rmdup)
    AAFTF taxpurge --phylum Ascomycota \
        -i $OUTDIR/${BASE}.vecscreen.fasta -o $OUTDIR/${BASE}.taxpurge.fasta

    AAFTF rmdup -c $CPU -i $OUTDIR/${BASE}.taxpurge.fasta -o $OUTDIR/${BASE}.rmdup.fasta

    AAFTF polish --method racon -c $CPU -lr $LONGREADS \
        -i $OUTDIR/${BASE}.rmdup.fasta -o $OUTDIR/${BASE}.racon.fasta

    AAFTF sort -i $OUTDIR/${BASE}.racon.fasta -o $OUTDIR/${BASE}.final.fasta

    AAFTF assess -i $OUTDIR/${BASE}.final.fasta -r $OUTDIR/${BASE}.stats.txt

.. _workflow-hybrid:

Illumina + long reads (pipeline_hybrid)
========================================

.. mermaid::

    flowchart TD
        illumina[/"Illumina FASTQ<br/>R1/R2 or single-end"/]
        longreads[/"long-read FASTQ<br/>ONT or PacBio"/]
        trim["<b>trim</b>"]
        filter["<b>filter</b>"]
        assemble["<b>assemble</b><br/>flye (long reads, default) or<br/>unicycler (Illumina + long reads)"]
        vecscreen["<b>vecscreen</b><br/>BLASTN vector + Euk/Prok/Mito screen"]
        taxpurge["<b>taxpurge</b> (with --phylum)<br/>sourmash taxonomy screen"]:::opt
        covpurge["<b>covpurge</b><br/>drop low-coverage contigs"]
        rmdup["<b>rmdup</b><br/>remove redundant contigs"]
        racon["<b>polish</b> (flye only)<br/>racon with the long reads"]
        short["<b>polish</b> (flye only)<br/>polypolish (paired) or<br/>pypolca (single-end)"]
        sort["<b>sort</b><br/>longest first, rename headers"]
        assess["<b>assess</b><br/>N50 / GC / telomere / gap stats"]

        illumina --> trim --> filter
        longreads --> assemble
        filter -. "unicycler" .-> assemble
        assemble --> vecscreen --> taxpurge --> covpurge --> rmdup --> racon --> short --> sort --> assess
        vecscreen -. "no --phylum" .-> covpurge
        rmdup -. "unicycler" .-> sort
        filter -. "filtered reads" .-> covpurge
        filter -. "filtered reads" .-> short
    classDef opt stroke-dasharray: 5 5

With ``--method unicycler``, Unicycler assembles both read types and polishes its own assembly, so
the polish steps are skipped.

.. list-table:: Files written (``{base}`` is the ``-o`` prefix)
   :header-rows: 1
   :widths: 14 30 30

   * - Step
     - Input
     - Output
   * - trim, filter
     - raw Illumina FASTQ
     - as in :ref:`workflow-short`
   * - assemble
     - long reads (flye), or filtered Illumina + long reads (unicycler)
     - ``{base}.flye.fasta`` or ``{base}.unicycler.fasta``
   * - vecscreen
     - assembled FASTA
     - ``{base}.vecscreen.fasta`` (+ ``{base}.mitochondria.fasta``)
   * - taxpurge (with ``--phylum``)
     - vecscreen FASTA
     - ``{base}.taxpurge.fasta``
   * - covpurge
     - taxpurge (or vecscreen) FASTA + filtered Illumina reads
     - ``{base}.covpurge.fasta``
   * - rmdup
     - covpurge FASTA
     - ``{base}.rmdup.fasta``
   * - polish (racon; flye only)
     - rmdup FASTA + long reads
     - ``{base}.racon.fasta``
   * - polish (polypolish / pypolca; flye only)
     - racon FASTA + filtered Illumina reads
     - ``{base}.polypolish.fasta`` (paired) or ``{base}.pypolca.fasta`` (single-end)
   * - sort
     - last polished FASTA (rmdup FASTA with unicycler)
     - ``{base}.final.fasta``
   * - assess
     - sorted FASTA
     - printed stats (optional ``--report`` file)

One command:

.. code-block:: bash

    AAFTF pipeline_hybrid \
        -1 reads/STRAINX_R1.fq.gz -2 reads/STRAINX_R2.fq.gz \
        -lr reads/STRAINX_ont.fq.gz \
        -o STRAINX -c 24 -m 96 --genome_size 40m \
        --phylum Ascomycota

Step by step (trim and filter the Illumina reads as in :ref:`workflow-short` first):

.. code-block:: bash

    CPU=24 BASE=STRAINX LONGREADS=reads/STRAINX_ont.fq.gz
    READ1=reads_trimmed/${BASE}_filtered_1.fastq.gz READ2=reads_trimmed/${BASE}_filtered_2.fastq.gz
    OUTDIR=genomes
    mkdir -p "$OUTDIR"

    AAFTF assemble --method flye -c $CPU \
        --longreads $LONGREADS --longread_type nano-hq --genome_size 40m \
        -o $OUTDIR/${BASE}.flye.fasta -w working_AAFTF/flye_${BASE}

    AAFTF vecscreen -c $CPU -i $OUTDIR/${BASE}.flye.fasta -o $OUTDIR/${BASE}.vecscreen.fasta

    # skip taxpurge if the organism's phylum is unknown (and pass the vecscreen FASTA to covpurge)
    AAFTF taxpurge --phylum Ascomycota \
        -i $OUTDIR/${BASE}.vecscreen.fasta -o $OUTDIR/${BASE}.taxpurge.fasta

    AAFTF covpurge -c $CPU --read1 $READ1 --read2 $READ2 \
        -i $OUTDIR/${BASE}.taxpurge.fasta -o $OUTDIR/${BASE}.covpurge.fasta

    AAFTF rmdup -c $CPU -i $OUTDIR/${BASE}.covpurge.fasta -o $OUTDIR/${BASE}.rmdup.fasta

    AAFTF polish --method racon -c $CPU -lr $LONGREADS \
        -i $OUTDIR/${BASE}.rmdup.fasta -o $OUTDIR/${BASE}.racon.fasta

    # polypolish needs paired reads; use --method pypolca with single-end reads
    AAFTF polish --method polypolish -c $CPU --read1 $READ1 --read2 $READ2 \
        -i $OUTDIR/${BASE}.racon.fasta -o $OUTDIR/${BASE}.polypolish.fasta

    AAFTF sort -i $OUTDIR/${BASE}.polypolish.fasta -o $OUTDIR/${BASE}.final.fasta

    AAFTF assess -i $OUTDIR/${BASE}.final.fasta -r $OUTDIR/${BASE}.stats.txt

Optional steps
==============

None of the ``pipeline_*`` subcommands run these; run them separately when needed:

* :doc:`commands/mito` -- assemble the mitochondrial genome from paired Illumina reads, then pass it
  to ``filter --screen_local`` to screen its reads out.
* :doc:`commands/fcs_screen` / :doc:`commands/fcs_gx_purge` -- NCBI FCS-based adaptor and
  cross-contamination screening, in place of or alongside :doc:`commands/taxpurge`.
* :doc:`commands/polish` -- other polishers (pypolca, POLCA, NextPolish2) or extra rounds.
* :doc:`commands/depth` -- per-contig read-depth report of the final assembly, with Illumina
  and/or long reads.

Each individual command page under :doc:`commands/index` documents its own defaults and cutoffs
in detail.
