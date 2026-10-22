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
(adaptor/vector BLAST screen) or, in its place, :doc:`commands/fcs_screen` (NCBI FCS-adaptor; both
screen for the same adaptors and vectors, so choose one), :doc:`commands/taxpurge` (sourmash taxonomy screen, run only
with ``-p/--phylum``; leave it out if the organism's taxonomy is unknown), :doc:`commands/covpurge`
(low Illumina read coverage; not in the long-read workflow) and :doc:`commands/rmdup` (redundant
contigs). Polishing, when the workflow has it, comes after cleaning so discarded contigs are never
polished.

In the diagrams below, the main column is what the ``pipeline_*`` subcommand runs by default;
dashed boxes and arrows are optional steps or alternative paths (e.g. taxpurge runs only with
``--phylum``).

.. _workflow-short:

Short reads (pipeline_short)
=============================

.. graphviz::
   :align: center

   digraph pipeline_short {
       graph [compound=true, bgcolor="white", pad="0.2", fontname="Helvetica", fontsize=11, nodesep=0.35, ranksep=0.32, newrank=true, labeljust=l, color="#9aa5b1", style="rounded,dashed"]
       node [shape=box, style="rounded,filled", fillcolor="#e3ecf9", color="#3d6aa8", fontname="Helvetica", fontsize=11, width=2.3, height=0.5, margin="0.12,0.05"]
       edge [color="#3d4852", arrowsize=0.7, fontname="Helvetica", fontsize=9, fontcolor="#5f6b7a"]
       reads [label="Illumina FASTQ\nR1/R2 or single-end", shape=note, fillcolor="#fdf1d8", color="#b7862d"]
       subgraph cluster_reads { label="Read processing"
           trim [label=<<b>trim</b><br/><font point-size="9">adapter + quality trimming</font>>]
           filter [label=<<b>filter</b><br/><font point-size="9">remove PhiX / vector / user contaminants</font>>]
       }
       mito [label=<<b>mito</b><br/><font point-size="9">optional, paired reads</font>>, style="rounded,dashed", fillcolor=white]
       assemble [label=<<b>assemble</b><br/><font point-size="9">spades (default) / megahit / unicycler</font>>]
       subgraph cluster_clean { label="Contig cleaning"
           vecscreen [label=<<b>vecscreen</b><br/><font point-size="9">BLAST adaptor/vector screen</font>>]
           covpurge [label=<<b>covpurge</b><br/><font point-size="9">drop low-coverage contigs</font>>]
           rmdup [label=<<b>rmdup</b><br/><font point-size="9">remove redundant contigs</font>>]
       }
       subgraph cluster_taxonomy { label="Taxonomic screen (optional)"
           taxpurge [label=<<b>taxpurge</b><br/><font point-size="9">sourmash, with -p</font>>, style="rounded,dashed", fillcolor=white]
       }
       fcs_screen [label=<<b>fcs_screen</b><br/><font point-size="9">NCBI FCS-adaptor, --fcs_screen<br/>(instead of vecscreen)</font>>, style="rounded,dashed", fillcolor=white]
       subgraph cluster_finish { label="Finishing"
           sort [label=<<b>sort</b><br/><font point-size="9">longest first, rename headers</font>>]
           depth [label=<<b>depth</b><br/><font point-size="9">coverage report, tag mito contigs</font>>]
           assess [label=<<b>assess</b><br/><font point-size="9">N50, GC, telomeres, gaps</font>>]
       }

       reads -> trim -> filter -> assemble -> vecscreen -> covpurge -> rmdup -> sort -> depth -> assess [weight=10]
       trim -> mito [style=dashed]
       mito -> filter [style=dashed, label=" --screen_local"]
       assemble -> fcs_screen [style=dashed]
       fcs_screen -> covpurge [style=dashed]
       {rank=same; vecscreen; fcs_screen}
       vecscreen -> taxpurge [style=dashed, lhead=cluster_taxonomy]
       fcs_screen -> taxpurge [style=dashed, lhead=cluster_taxonomy, constraint=false]
       taxpurge -> covpurge [style=dashed, ltail=cluster_taxonomy]
       {rank=same; filter; mito}
       {rank=same; covpurge; taxpurge}
   }

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
   * - vecscreen (or fcs_screen with ``--fcs_screen``)
     - assembled FASTA
     - ``{base}.vecscreen.fasta`` (or ``{base}.fcs_screen.fasta``)
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
     - ``{base}.sorted.fasta``
   * - depth
     - sorted FASTA + reads
     - ``{base}.final.fasta`` (likely mitochondrial contigs tagged ``[location=mitochondrion]``),
       ``{base}.coverage_stats.txt``, ``{base}.coverage_stats.mito_candidates.tsv`` (+ plots)
   * - assess
     - final FASTA
     - printed stats (optional ``--report`` file)

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
    # or, instead of vecscreen (not both): AAFTF fcs_screen -i ... -o $OUTDIR/${BASE}.fcs_screen.fasta

    # skip taxpurge if the organism's phylum is unknown (and pass the vecscreen FASTA to covpurge)
    AAFTF taxpurge --phylum Ascomycota \
        -i $OUTDIR/${BASE}.vecscreen.fasta -o $OUTDIR/${BASE}.taxpurge.fasta

    AAFTF covpurge -c $CPU \
        -i $OUTDIR/${BASE}.taxpurge.fasta -o $OUTDIR/${BASE}.covpurge.fasta \
        --read1 $TRIMREAD/${BASE}_filtered_1.fastq.gz --read2 $TRIMREAD/${BASE}_filtered_2.fastq.gz

    AAFTF rmdup -c $CPU -i $OUTDIR/${BASE}.covpurge.fasta -o $OUTDIR/${BASE}.rmdup.fasta

    AAFTF sort -i $OUTDIR/${BASE}.rmdup.fasta -o $OUTDIR/${BASE}.sorted.fasta

    # coverage report; --tagged writes the final assembly with likely mitochondrial contigs tagged
    AAFTF depth -i $OUTDIR/${BASE}.sorted.fasta \
        --read1 $TRIMREAD/${BASE}_filtered_1.fastq.gz --read2 $TRIMREAD/${BASE}_filtered_2.fastq.gz \
        -c $CPU -o $OUTDIR/${BASE}.coverage_stats.txt --tagged $OUTDIR/${BASE}.final.fasta

    AAFTF assess -i $OUTDIR/${BASE}.final.fasta -r $OUTDIR/${BASE}.stats.txt

.. _workflow-long:

Long reads (pipeline_long)
===========================

.. graphviz::
   :align: center

   digraph pipeline_long {
       graph [compound=true, bgcolor="white", pad="0.2", fontname="Helvetica", fontsize=11, nodesep=0.35, ranksep=0.32, newrank=true, labeljust=l, color="#9aa5b1", style="rounded,dashed"]
       node [shape=box, style="rounded,filled", fillcolor="#e3ecf9", color="#3d6aa8", fontname="Helvetica", fontsize=11, width=2.3, height=0.5, margin="0.12,0.05"]
       edge [color="#3d4852", arrowsize=0.7, fontname="Helvetica", fontsize=9, fontcolor="#5f6b7a"]
       reads [label="long-read FASTQ\nONT or PacBio", shape=note, fillcolor="#fdf1d8", color="#b7862d"]
       assemble [label=<<b>assemble</b><br/><font point-size="9">flye</font>>]
       subgraph cluster_clean { label="Contig cleaning"
           vecscreen [label=<<b>vecscreen</b><br/><font point-size="9">BLAST adaptor/vector screen</font>>]
           rmdup [label=<<b>rmdup</b><br/><font point-size="9">remove redundant contigs</font>>]
       }
       subgraph cluster_taxonomy { label="Taxonomic screen (optional)"
           taxpurge [label=<<b>taxpurge</b><br/><font point-size="9">sourmash, with -p</font>>, style="rounded,dashed", fillcolor=white]
       }
       fcs_screen [label=<<b>fcs_screen</b><br/><font point-size="9">NCBI FCS-adaptor, --fcs_screen<br/>(instead of vecscreen)</font>>, style="rounded,dashed", fillcolor=white]
       subgraph cluster_polish { label="Polishing"
           racon [label=<<b>polish</b><br/><font point-size="9">racon with the long reads</font>>]
       }
       subgraph cluster_finish { label="Finishing"
           sort [label=<<b>sort</b><br/><font point-size="9">longest first, rename headers</font>>]
           depth [label=<<b>depth</b><br/><font point-size="9">coverage report, tag mito contigs</font>>]
           assess [label=<<b>assess</b><br/><font point-size="9">N50, GC, telomeres, gaps</font>>]
       }

       reads -> assemble -> vecscreen -> rmdup -> racon -> sort -> depth -> assess [weight=10]
       assemble -> fcs_screen [style=dashed]
       fcs_screen -> rmdup [style=dashed]
       {rank=same; vecscreen; fcs_screen}
       vecscreen -> taxpurge [style=dashed, lhead=cluster_taxonomy]
       fcs_screen -> taxpurge [style=dashed, lhead=cluster_taxonomy, constraint=false]
       taxpurge -> rmdup [style=dashed, ltail=cluster_taxonomy]
       {rank=same; rmdup; taxpurge}
   }

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
   * - vecscreen (or fcs_screen with ``--fcs_screen``)
     - flye FASTA
     - ``{base}.vecscreen.fasta`` (or ``{base}.fcs_screen.fasta``)
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
     - ``{base}.sorted.fasta``
   * - depth
     - sorted FASTA + reads
     - ``{base}.final.fasta`` (likely mitochondrial contigs tagged ``[location=mitochondrion]``),
       ``{base}.coverage_stats.txt``, ``{base}.coverage_stats.mito_candidates.tsv`` (+ plots)
   * - assess
     - final FASTA
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

    AAFTF sort -i $OUTDIR/${BASE}.racon.fasta -o $OUTDIR/${BASE}.sorted.fasta

    AAFTF depth -i $OUTDIR/${BASE}.sorted.fasta -lr $LONGREADS \
        -c $CPU -o $OUTDIR/${BASE}.coverage_stats.txt --tagged $OUTDIR/${BASE}.final.fasta

    AAFTF assess -i $OUTDIR/${BASE}.final.fasta -r $OUTDIR/${BASE}.stats.txt

.. _workflow-hybrid:

Illumina + long reads (pipeline_hybrid)
========================================

.. graphviz::
   :align: center

   digraph pipeline_hybrid {
       graph [compound=true, bgcolor="white", pad="0.2", fontname="Helvetica", fontsize=11, nodesep=0.35, ranksep=0.32, newrank=true, labeljust=l, color="#9aa5b1", style="rounded,dashed"]
       node [shape=box, style="rounded,filled", fillcolor="#e3ecf9", color="#3d6aa8", fontname="Helvetica", fontsize=11, width=2.3, height=0.5, margin="0.12,0.05"]
       edge [color="#3d4852", arrowsize=0.7, fontname="Helvetica", fontsize=9, fontcolor="#5f6b7a"]
       illumina [label="Illumina FASTQ\nR1/R2 or single-end", shape=note, fillcolor="#fdf1d8", color="#b7862d"]
       longreads [label="long-read FASTQ\nONT or PacBio", shape=note, fillcolor="#fdf1d8", color="#b7862d"]
       subgraph cluster_reads { label="Read processing"
           trim [label=<<b>trim</b><br/><font point-size="9">adapter + quality trimming</font>>]
           filter [label=<<b>filter</b><br/><font point-size="9">remove PhiX / vector / user contaminants</font>>]
       }
       assemble [label=<<b>assemble</b><br/><font point-size="9">flye (default): long reads<br/>unicycler: both read types</font>>]
       subgraph cluster_clean { label="Contig cleaning"
           vecscreen [label=<<b>vecscreen</b><br/><font point-size="9">BLAST adaptor/vector screen</font>>]
           covpurge [label=<<b>covpurge</b><br/><font point-size="9">drop low-coverage contigs</font>>]
           rmdup [label=<<b>rmdup</b><br/><font point-size="9">remove redundant contigs</font>>]
       }
       subgraph cluster_taxonomy { label="Taxonomic screen (optional)"
           taxpurge [label=<<b>taxpurge</b><br/><font point-size="9">sourmash, with -p</font>>, style="rounded,dashed", fillcolor=white]
       }
       fcs_screen [label=<<b>fcs_screen</b><br/><font point-size="9">NCBI FCS-adaptor, --fcs_screen<br/>(instead of vecscreen)</font>>, style="rounded,dashed", fillcolor=white]
       subgraph cluster_polish { label="Polishing (flye only)"
           racon [label=<<b>polish</b><br/><font point-size="9">racon with the long reads</font>>]
           short [label=<<b>polish</b><br/><font point-size="9">polypolish (paired) / pypolca (single-end)</font>>]
       }
       subgraph cluster_finish { label="Finishing"
           sort [label=<<b>sort</b><br/><font point-size="9">longest first, rename headers</font>>]
           depth [label=<<b>depth</b><br/><font point-size="9">coverage report, tag mito contigs</font>>]
           assess [label=<<b>assess</b><br/><font point-size="9">N50, GC, telomeres, gaps</font>>]
       }

       illumina -> trim -> filter [weight=10]
       longreads -> assemble [weight=10]
       filter -> assemble [style=dashed, label=" unicycler"]
       assemble -> vecscreen -> covpurge -> rmdup -> racon -> short -> sort -> depth -> assess [weight=10]
       assemble -> fcs_screen [style=dashed]
       fcs_screen -> covpurge [style=dashed]
       {rank=same; vecscreen; fcs_screen}
       vecscreen -> taxpurge [style=dashed, lhead=cluster_taxonomy]
       fcs_screen -> taxpurge [style=dashed, lhead=cluster_taxonomy, constraint=false]
       taxpurge -> covpurge [style=dashed, ltail=cluster_taxonomy]
       rmdup -> sort [style=dashed, label=" unicycler", constraint=false]
       {rank=same; illumina; longreads}
       {rank=same; covpurge; taxpurge}
   }

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
   * - vecscreen (or fcs_screen with ``--fcs_screen``)
     - assembled FASTA
     - ``{base}.vecscreen.fasta`` (or ``{base}.fcs_screen.fasta``)
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
     - ``{base}.sorted.fasta``
   * - depth
     - sorted FASTA + reads
     - ``{base}.final.fasta`` (likely mitochondrial contigs tagged ``[location=mitochondrion]``),
       ``{base}.coverage_stats.txt``, ``{base}.coverage_stats.mito_candidates.tsv`` (+ plots)
   * - assess
     - final FASTA
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

    AAFTF sort -i $OUTDIR/${BASE}.polypolish.fasta -o $OUTDIR/${BASE}.sorted.fasta

    AAFTF depth -i $OUTDIR/${BASE}.sorted.fasta --read1 $READ1 --read2 $READ2 \
        -c $CPU -o $OUTDIR/${BASE}.coverage_stats.txt --tagged $OUTDIR/${BASE}.final.fasta

    AAFTF assess -i $OUTDIR/${BASE}.final.fasta -r $OUTDIR/${BASE}.stats.txt

Optional steps
==============

None of the ``pipeline_*`` subcommands run these by default:

* :doc:`commands/mito` -- assemble the mitochondrial genome from paired Illumina reads, then pass it
  to ``filter --screen_local`` to screen its reads out.
* :doc:`commands/fcs_screen` -- NCBI FCS-adaptor adaptor/vector screen; the pipelines run it
  instead of vecscreen with ``--fcs_screen``.
* :doc:`commands/fcs_gx_purge` -- NCBI FCS-GX cross-contamination purge. It needs its whole
  database (~470 GB) in RAM, so it is not a pipeline step: run it separately, on a large-memory
  node, on one or more finished assemblies.
* :doc:`commands/polish` -- other polishers (pypolca, POLCA, NextPolish2) or extra rounds.
* :doc:`commands/depth` -- the pipelines already run it to tag mitochondrial contigs; run it again
  with other reads (e.g. both Illumina and long reads) or settings for another coverage report.

Each individual command page under :doc:`commands/index` documents its own defaults and cutoffs
in detail.
