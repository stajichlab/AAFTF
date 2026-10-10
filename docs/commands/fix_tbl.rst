=======
fix_tbl
=======

Adjusts feature coordinates in an NCBI ``.tbl`` annotation table after :doc:`fcs_screen` (or any
NCBI FCS trimming step) has trimmed bases off the start/end of one or more sequences -- otherwise
annotated feature coordinates would point past the end of, or into the wrong part of, the trimmed
sequence.

Algorithm
=========

1. Parses the input ``.tbl`` file into per-sequence features (``>Feature seqid`` blocks; each
   feature is a ``start<TAB>end<TAB>key`` line, optional extra ``start<TAB>end`` interval lines,
   and its qualifier lines). Minus-strand features (``start > end``) and ``<``/``>`` partial
   markers are kept as written.
2. Parses the NCBI FCS action report (tab-separated: accession, length, action, range(s), ...),
   reading every row after the ``#accession`` header line.
3. For each sequence, works out which original bases remain: trims starting at position 1 remove
   bases from the **left**, trims ending at the original length remove bases from the **right**
   (the furthest trim wins when there are several). A sequence marked ``ACTION_EXCLUDE``, or
   trimmed completely, is dropped from the table with a warning. Trims inside a contig cannot be
   fixed by shifting coordinates and are reported as a warning.
4. Clips every feature to the remaining bases and shifts it so the first remaining base becomes 1.
   A feature cut at its 5' end gets ``<`` on its first coordinate, and one cut at its 3' end gets
   ``>`` on its last coordinate (on either strand, including when a whole interval of a
   multi-interval feature was trimmed away). Features that lie entirely inside a trimmed region
   are dropped with a warning.
5. Writes the corrected ``.tbl`` file to ``{output}.tmp`` and renames it to ``-o/--output`` only
   once complete, so a malformed input never leaves a half-written output. Malformed ``.tbl``
   lines (a header other than ``>Feature``, a feature before the first header) raise
   ``ValueError`` (exit code 1).

Invocation
==========

.. code-block:: text

    AAFTF fix_tbl -t TABLE -r FILE -o OUTPUT [-q] [-v]

All three of ``-t/--table`` (alias ``--infile``; original ``.tbl``), ``-r/--report`` (NCBI FCS
action report) and ``-o/--output`` (corrected ``.tbl``) are required, and all are plain file
paths (``-`` is not treated as stdin/stdout). ``./fix_tbl.log`` is written only with ``-v``.

Example
=======

.. code-block:: bash

    AAFTF fix_tbl -t STRAINX.tbl -r STRAINX.fcs_adaptor_report.txt -o STRAINX.fixed.tbl

Typical use: after running :doc:`fcs_screen` on a contig set that already has NCBI-format
annotation (``.tbl``) generated against the *untrimmed* sequences, run ``fix_tbl`` to bring the
annotation coordinates back in sync with the trimmed FASTA before resubmission.
