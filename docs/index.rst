======================================================================================================
**AAFTF**: *Automatic Assembly For The Fungi. a utility for Genome Assembly, Cleanup, and Assessment.*
======================================================================================================

AAFTF is a menu-driven command line toolkit (``AAFTF <subcommand> [options]``) that takes raw
paired-end Illumina reads (optionally with PacBio/ONT long reads) through quality trimming,
contamination filtering, assembly, vector/contaminant screening, duplicate removal, polishing,
sorting, and quality assessment -- producing a clean, submission-ready fungal genome assembly plus
summary statistics describing assembly completeness and sequencing success.

Every subcommand is invoked the same way::

    AAFTF <subcommand> [options]

and every subcommand accepts at minimum ``-v/--debug`` (verbose logging, keep temp files) and
``--pipe`` (suppress "next command" hints; used internally when ``AAFTF pipeline`` chains steps
together). Many subcommands have one or more shorter aliases (e.g. ``AAFTF asm`` for
``AAFTF assemble``, ``AAFTF stats`` for ``AAFTF assess``) -- see each command page for its aliases.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   workflow
   commands/index

Quick links
===========

* :doc:`installation` -- installing via conda/pixi, or running from the Docker/Singularity image
* :doc:`workflow` -- the recommended raw-reads-to-assembly pipeline and what each step produces
* :doc:`commands/index` -- reference documentation for every subcommand, with algorithm details,
  cutoffs/defaults, and worked examples

Source, issues, and citation
=============================

* GitHub: https://github.com/stajichlab/AAFTF
* Citation: Palmer JM and Stajich JE. (2023). *Automatic assembly for the fungi (AAFTF): genome
  assembly pipeline* (v0.5.0). Zenodo. doi: 10.5281/zenodo.1620526
