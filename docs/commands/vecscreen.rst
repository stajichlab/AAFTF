=========
vecscreen
=========

Screens assembled contigs for cloning-vector/adapter contamination (replicating NCBI's VecScreen
service) and for common eukaryotic/prokaryotic/mitochondrial contaminant sequences that slipped
through the read-level :doc:`filter` step, using BLASTN.

Algorithm
=========

Runs three successive BLASTN screens, in this order, against BLAST databases built in the work
directory from the ``univec``, ``euks``, ``proks`` and ``mitodb`` databases installed by
``AAFTF database`` (the run stops with exit code 2 if any of them is missing):

1. **CONTAM_EUKS / CONTAM_PROKS screen.** BLASTN of the assembly against NCBI's
   ``contam_in_euks``/``contam_in_prok`` reference sets (``-dust yes -soft_masking true``). A hit
   region is trimmed out of its contig (splitting the contig into fragments around the hit) if it
   meets any of: **>=98% identity over >=50 bp**, **>=94% identity over >=100 bp**, or
   **>=90% identity over >=200 bp** (``-perc_identity`` cutoff itself defaults to 90.0, overridable
   with ``-pid/--percent_id``).
2. **Mitochondrial screen.** BLASTN of the (euk/prok-cleaned) assembly against the RefSeq
   mitochondrion database at **>=98.6% identity**; any hit **>=120 bp** flags the *entire contig*
   for removal from the nuclear assembly (written instead to a separate
   ``{prefix}.mitochondria.fasta`` file).
3. **VecScreen (UniVec) screen.** Iterative BLASTN against UniVec using NCBI's documented
   VecScreen scoring parameters (``-reward 1 -penalty -5 -gapopen 3 -gapextend 3 -evalue 700``),
   classifying each hit as terminal (within 25 bp of a contig end) or internal, and as
   weak/moderate/strong per NCBI's VecScreen score thresholds:

   .. list-table::
      :header-rows: 1
      :widths: 25 25 25 25

      * - Match strength
        - Terminal score
        - Internal score
        - Random-match rate
      * - Strong
        - >= 24
        - >= 30
        - 1 in 1,000,000 (350 kb query)
      * - Moderate
        - 19-23
        - 25-29
        - 1 in 1,000
      * - Weak
        - 16-18
        - 23-24
        - 1 in 40

   Weak matches are always ignored. ``-s/--stringency high`` (default) acts on moderate+strong
   matches; ``low`` acts on strong matches only. Terminal hits are trimmed off the nearest end; internal hits split the
   contig into pieces around the hit. This repeats in rounds against the progressively-cleaned
   assembly until no further vector hits are found. Fragments shorter than 200 bp after trimming
   are dropped.

Contigs removed entirely by the mito screen are written to ``{prefix}.mitochondria.fasta``
(``{prefix}`` = the output name up to its last ``.f...`` extension), in the same directory as
``-o/--outfile``; the cleaned assembly is written to ``-o/--outfile`` itself. BLAST
reports and intermediate FASTA files go to the work directory, together with
``vecscreen.log`` (the work directory is kept when given with ``-w``, with ``-v``, or when the
run fails).

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-pid/--percent_id``
     - 90.0
     - BLASTN ``-perc_identity`` for the Euk/Prok contamination screen
   * - Mito match identity
     - 98.6% (hard-coded)
     - BLASTN identity cutoff for the mitochondrial screen
   * - Mito match length
     - 120 bp (hard-coded)
     - Minimum alignment length to flag a whole contig as mitochondrial
   * - ``-s/--stringency``
     - high
     - VecScreen sensitivity: ``high`` (moderate+strong) or ``low`` (strong only)
   * - Minimum surviving fragment
     - 200 bp (hard-coded)
     - Trimmed/split fragments shorter than this are dropped

Invocation
==========

.. code-block:: text

    AAFTF vecscreen -i FASTA -o FASTA [-w DIR] [-pid PCT] [-s {high,low}]
                    [-c INT] [-q] [-v]

``-i/--input`` (alias ``--infile``; assembly FASTA) and ``-o/--outfile`` are required. Unless
``-q`` is given, the suggested next command (:doc:`taxpurge`) is printed at the end.

Example
=======

.. code-block:: bash

    # writes genomes/STRAINX.vecscreen.fasta and genomes/STRAINX.vecscreen.mitochondria.fasta
    AAFTF vecscreen -c 16 -i genomes/STRAINX.spades.fasta -o genomes/STRAINX.vecscreen.fasta

Next step: :doc:`taxpurge` (sourmash-based) or :doc:`fcs_gx_purge` (NCBI FCS-GX-based) --
alternative approaches to the same taxonomic-contamination-purge step.
