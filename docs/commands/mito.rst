====
mito
====

De novo assembles the mitochondrial genome from trimmed paired-end Illumina reads, independently
of the nuclear genome assembly, using NOVOPlasty's seed-and-extend organelle assembler.

Algorithm
=========

1. Subsamples the reads (``--subsample``, default 1,500,000 pairs; ``0`` keeps all reads): BBTools
   ``reformat.sh`` keeps that many randomly chosen read pairs, with a fixed seed so reruns keep the
   same pairs. Mitochondrial reads are usually at far higher coverage than nuclear ones, so this
   saves time and memory without starving the assembly. (NOVOPlasty may also subsample on its
   own to stay within ``-m``.)
   Then estimates read length from the FASTQ (``estimate_read_length``).
2. Selects the **seed**, the sequence NOVOPlasty starts assembling from and extends with matching
   reads: a user-supplied ``--seed`` FASTA (e.g. a gene, or the mitochondrial genome of a
   related species) or, by default, a bundled *Aspergillus nidulans* cytochrome-b (cob) fragment
   (``aaftf/data/mito-seed.fasta``).
3. Writes a NOVOPlasty config file from the ``aaftf/data/novoplasty-config.txt`` template
   (substituting project name, min/max genome length, max memory, seed, read length, and the
   forward/reverse FASTQ paths) and runs ``NOVOPlasty.pl -c novo-config.txt``.
4. Parses NOVOPlasty's output directory for (in priority order) a
   ``Circularized_assembly_*``, ``Contigs_1_*``, or ``Uncircularized_assemblies_*`` file.
5. If circularized, rotates/reorients the genome to start at the cytochrome-b (cob) gene: aligns
   a bundled spoa consensus of fungal cob genes (``aaftf/data/mito-start-cob.fasta``) against the
   assembly with
   ``minimap2 -x map-ont``, and rotates the sequence to begin at that alignment's start
   coordinate (reverse-complementing if the hit is on the minus strand). Rotation is skipped
   (contig kept as-is, with a warning) if zero or multiple alignments are found, or if the
   computed rotation offset falls outside the sequence length.
6. If not circularized, all contigs are renumbered and concatenated into the output FASTA as-is,
   with a warning that circularization failed.

Requires ``NOVOPlasty.pl`` and ``minimap2`` (plus ``reformat.sh`` when subsampling) on ``$PATH``;
a missing tool raises an error (exit code 2) before any work starts. If NOVOPlasty produces no
assembly, a ``RuntimeError`` is raised (exit code 1).

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``--minlen``
     - 10000
     - Minimum expected mitochondrial genome size (NOVOPlasty search bound)
   * - ``--maxlen``
     - 100000
     - Maximum expected mitochondrial genome size (NOVOPlasty search bound)
   * - ``--subsample``
     - 1500000
     - Number of read pairs to keep before assembling; ``0`` keeps all reads
   * - ``-m/--memory``
     - 8 (GB)
     - NOVOPlasty max RAM, passed into the generated config
   * - ``--seed``
     - bundled *A. nidulans* cob fragment (``aaftf/data/mito-seed.fasta``)
     - Where NOVOPlasty starts assembling (the sequence it extends)

Invocation
==========

.. code-block:: text

    AAFTF mito -1 FASTQ -2 FASTQ [-o FASTA] [--minlen BP] [--maxlen BP]
               [-s FASTA] [--subsample PAIRS] [-m GB]
               [-w DIR] [-q] [-v]

``-1/--read1`` and ``-2/--read2`` are required (``mito`` only supports paired-end data); ``-o/--out``
defaults to ``mito.fasta``.

Example
=======

.. code-block:: bash

    AAFTF mito -1 reads_trimmed/STRAINX_1P.fastq.gz -2 reads_trimmed/STRAINX_2P.fastq.gz \
        -o STRAINX.mito.fasta

    # Seed with a related species' mitogenome
    AAFTF mito -1 STRAINX_1P.fastq.gz -2 STRAINX_2P.fastq.gz -o STRAINX.mito.fasta \
        --seed related_species_mito.fasta

``mito`` is not part of ``pipeline_short``, ``pipeline_long`` or ``pipeline_hybrid``. To keep mitochondrial reads out of the nuclear assembly,
run it before :doc:`filter` and pass its output with ``-s``; unless ``-q`` is given, ``mito``
prints this next command, e.g.
``AAFTF filter -1 STRAINX_1P.fastq.gz -2 STRAINX_2P.fastq.gz -s STRAINX.mito.fasta -o STRAINX``.
