=======
fix_tbl
=======

**Aliases:** ``fix``

Adjusts feature coordinates in an NCBI ``.tbl`` annotation table after :doc:`fcs_screen` (or any
NCBI FCS trimming step) has trimmed bases off the start/end of one or more sequences -- otherwise
annotated feature coordinates would point past the end of, or into the wrong part of, the trimmed
sequence.

Algorithm
=========

1. Parses the input ``.tbl`` file into per-sequence feature lists (``>Feature seqid`` blocks,
   each containing start/end/feature-type rows and qualifier lines).
2. Parses the NCBI FCS trim-adjustment report (5-column TSV: accession, original length, action,
   trimmed range(s), ...), keeping only ``ACTION_TRIM`` rows.
3. For each affected sequence, determines whether the trim was at the **left** end (trim range
   starts at position 1) or the **right** end (trim range ends at the original sequence length).
   Trims spanning the interior of a contig are not auto-correctable and are reported as a warning
   instead of silently mis-adjusting coordinates.
4. Shifts every feature's start/end coordinates by the trimmed amount (clamping to a minimum of 1
   for a left trim; clamping the end to the new right boundary for a right trim), preserving
   NCBI's partial-feature markers (``<``/``>`` prefixes) on the coordinates.
5. Writes the corrected ``.tbl`` file.

Invocation
==========

.. code-block:: text

    AAFTF fix_tbl -t TABLE -r REPORT -o OUTPUT [-v] [--pipe]

All three of ``-t/--table`` (original ``.tbl``), ``-r/--report`` (NCBI FCS adjustment report),
and ``-o/--output`` (corrected ``.tbl``) are required.

Example
=======

.. code-block:: bash

    AAFTF fix_tbl -t STRAINX.tbl -r STRAINX.fcs_adaptor_report.txt -o STRAINX.fixed.tbl

Typical use: after running :doc:`fcs_screen` on a contig set that already has NCBI-format
annotation (``.tbl``) generated against the *untrimmed* sequences, run ``fix_tbl`` to bring the
annotation coordinates back in sync with the trimmed FASTA before resubmission.
