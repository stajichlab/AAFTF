=========
vecscreen
=========

**Aliases:** ``vectorscreen``, ``vector_blast``

Screens assembled contigs for cloning-vector/adapter contamination (replicating NCBI's VecScreen
service) and for common eukaryotic/prokaryotic/mitochondrial contaminant sequences that slipped
through the read-level :doc:`filter` step, using BLASTN.

Algorithm
=========

Runs three successive BLASTN screens, in this order, against BLAST databases built from
``resources.DB_Links`` (downloaded/cached under ``$AAFTF_DB`` or the working directory):

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

   ``-s/--stringency high`` (default) keeps/acts on moderate+strong matches; ``low`` acts on
   strong matches only. Terminal hits are trimmed off the nearest end; internal hits split the
   contig into pieces around the hit. This repeats in rounds against the progressively-cleaned
   assembly until no further vector hits are found. Fragments shorter than 200 bp after trimming
   are dropped.

Contigs removed entirely by the mito screen are written to ``{prefix}.mitochondria.fasta``
alongside the main cleaned output.

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

    AAFTF vecscreen -i INFILE -o OUTFILE [-c CPUS] [-pid PERCENT_ID]
                    [-s {high,low}] [--AAFTF_DB DIR] [-w WORKDIR] [-v] [--pipe]

``-i/--input`` (assembly FASTA) and ``-o/--outfile`` are required.

Example
=======

.. code-block:: bash

    AAFTF vecscreen -c 16 -i genomes/STRAINX.spades.fasta -o genomes/STRAINX.vecscreen.fasta

Next step: :doc:`sourpurge` (sourmash-based) or :doc:`fcs_gx_purge` (NCBI FCS-GX-based) --
alternative approaches to the same taxonomic-contamination-purge step.
