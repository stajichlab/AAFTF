====
mito
====

**Aliases:** ``mito_asm``, ``mitochondria``

De novo assembles the mitochondrial genome from trimmed paired-end Illumina reads, independently
of the nuclear genome assembly, using NOVOPlasty's seed-and-extend organelle assembler.

Algorithm
=========

1. Estimates read length from the input FASTQ (``GuessRL``).
2. Selects a seed sequence: a user-supplied ``--seed`` FASTA, a ``--reference`` mitochondrial
   genome (also enables NOVOPlasty's reference-guided mode), or -- by default -- a bundled
   *Aspergillus nidulans* cytochrome-b (COB) fragment (``AAFTF/data/mito-seed.fasta``).
3. Writes a NOVOPlasty config file from the ``AAFTF/data/novoplasty-config.txt`` template
   (substituting project name, min/max genome length, max memory, seed, read length, and the
   forward/reverse FASTQ paths) and runs ``NOVOPlasty.pl -c novo-config.txt``.
4. Parses NOVOPlasty's output directory for (in priority order) a
   ``Circularized_assembly_*``, ``Contigs_1_*``, or ``Uncircularized_assemblies_*`` file.
5. If circularized, rotates/reorients the genome to start at a chosen gene: aligns a start
   sequence (``--starting``, default cytochrome b) against the assembly with
   ``minimap2 -x map-ont``, and rotates the sequence to begin at that alignment's start
   coordinate (reverse-complementing if the hit is on the minus strand). Rotation is skipped
   (contig kept as-is, with a warning) if zero or multiple alignments are found, or if the
   computed rotation offset falls outside the sequence length.
6. If not circularized, all contigs are renumbered and concatenated into the output FASTA as-is,
   with a warning that circularization failed.

Requires ``NOVOPlasty.pl`` and ``minimap2`` on ``$PATH``; exits immediately if either is missing.

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
   * - NOVOPlasty max RAM
     - 75% of detected system RAM
     - Passed into the generated config
   * - ``--seed``
     - bundled *A. nidulans* COB fragment
     - Seed read for NOVOPlasty's extension algorithm
   * - ``--starting``
     - same COB fragment
     - Sequence used to rotate/orient the final circular genome

Invocation
==========

.. code-block:: text

    AAFTF mito -l LEFT -r RIGHT -o OUT [--minlen N] [--maxlen N]
               [-s SEED] [--starting FASTA] [--reference REFGENOME]
               [-w WORKDIR] [-v] [--pipe]

``-l/--left``, ``-r/--right``, and ``-o/--out`` are required; ``mito`` only supports paired-end
data.

Example
=======

.. code-block:: bash

    AAFTF mito -l reads_trimmed/STRAINX_1P.fastq.gz -r reads_trimmed/STRAINX_2P.fastq.gz \
        -o STRAINX.mito.fasta

    # Reference-guided assembly against a related species' mitogenome
    AAFTF mito -l STRAINX_1P.fastq.gz -r STRAINX_2P.fastq.gz -o STRAINX.mito.fasta \
        --reference related_species_mito.fasta

In the full ``pipeline`` workflow, the resulting ``{base}.mito.fasta`` is fed to :doc:`filter` as
an extra ``--screen_local`` contamination reference, so mitochondrial reads assembled here are
kept separate from (and not double-counted against) the nuclear assembly.
