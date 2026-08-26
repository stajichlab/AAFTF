======
assess
======

**Aliases:** ``stats``

Computes standard genome assembly completeness/summary statistics -- the final QC report
describing what the pipeline produced. Pure BioPython; no external tools required.

Algorithm
=========

For each contig, tallies GC+S content, N-base count and N-gap count (runs of ``N`` characters,
counted via a regex), and soft-masked (lowercase) base count. From the sorted length distribution
it computes:

* **N50/L50** and **N90/L90** -- length (N) and count (L) of the contig at which cumulative
  assembly length (summed longest-to-shortest) first reaches 50%/90% of total assembly size.
* **MIN / MAX / MEDIAN / MEAN** contig length, **CONTIG COUNT**, **TOTAL LENGTH**, **GC%**.
* **N GAP COUNT** / **TOTAL N BASES** -- number of distinct N-runs and total N bases (assembly
  gaps).
* **BASES MASKED / PERCENT MASKED** -- only reported if the assembly contains any lowercase
  (soft-masked) bases.
* **Telomere detection** -- for each contig, scans the first and last ``--telomere_window`` bp
  (default 200) for >= ``-n/--telomere_n_repeat`` (default 2) tandem copies of the
  ``-t/--telomere_monomer`` repeat motif (default ``TAA[C]+``, the canonical fungal telomere
  repeat) at the 5' end, and its reverse complement at the 3' end. Reports counts of
  **TELOMERE FWD**, **TELOMERE REV**, and **T2T SCAFFOLDS** (contigs with a telomere repeat found
  at *both* ends -- fully telomere-to-telomere assembled chromosomes/scaffolds).

Output is printed to stdout and, if ``-r/--report`` is given, also written to that file.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-t/--telomere_monomer``
     - ``TAA[C]+``
     - Telomere repeat motif to search for (regex; default matches the fungal ``TTAGGG``-family
       repeat allowing extra ``C``\ s)
   * - ``-n/--telomere_n_repeat``
     - 2
     - Minimum tandem repeat count required to call a telomere
   * - ``--telomere_window``
     - 200 (bp)
     - Window scanned at each contig end for telomere repeats

Invocation
==========

.. code-block:: text

    AAFTF assess -i INPUT [-r REPORT] [-t TELOMERE_MONOMER] [-n TELOMERE_N_REPEAT]
                 [--telomere_window N] [-v] [--pipe]

``-i/--input`` is required.

Example
=======

.. code-block:: bash

    AAFTF assess -i genomes/STRAINX.final.fasta -r genomes/STRAINX.stats.txt

    # non-fungal telomere motif, wider search window
    AAFTF assess -i genomes/STRAINX.final.fasta -t 'TTAGGG' --telomere_window 500

Optional next step: :doc:`depth` for a per-contig read coverage report.
