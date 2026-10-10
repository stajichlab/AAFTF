==========
dependency
==========

Checks that the external tools and Python packages AAFTF uses are installed, and reports which
are missing and which subcommands need them. Run it after installing, or when a step fails because
a tool is not on ``$PATH``.

Algorithm
=========

Each tool is looked up on ``$PATH`` and each Python package is imported. They are checked in two
groups, and printed in the order listed below:

* **Needed by the default pipelines** (``REQUIRED_TOOLS`` / ``REQUIRED_PYTHON_PACKAGES`` in
  ``aaftf/dependency.py``) -- what ``pipeline_short``, ``pipeline_long`` and ``pipeline_hybrid``
  need with their default settings: ``bbduk.sh``, ``shuffle.sh``, ``reformat.sh``, ``java``,
  ``spades.py``, ``blastn``, ``makeblastdb``, ``bwa``, ``samtools``, ``minimap2``,
  ``flye``, ``racon``, ``polypolish``, ``pypolca``, ``freebayes``; Python packages
  ``biopython`` and ``psutil``.
* **Optional steps and non-default options** (``OPTIONAL_TOOLS`` / ``OPTIONAL_PYTHON_PACKAGES``):
  ``pigz``, ``sourmash`` (``taxpurge``), ``fastp``, ``trimmomatic``, ``bowtie2``, ``bowtie2-build``, ``megahit``,
  ``unicycler``, ``aaftf_polca.sh``, ``nextPolish2``, ``yak``, ``mosdepth``, ``NOVOPlasty.pl``, ``run_fcsadaptor.sh``, ``singularity``, ``apptainer``,
  ``docker``, ``run_gx.py``; Python packages ``matplotlib`` (``depth`` plots) and
  ``sourmash_plugin_branchwater`` (``taxpurge``).

For each tool an ``[OK]`` line shows where it was found, and a ``[MISSING]`` line shows which
subcommands/options use it (e.g. ``used by: filter (--aligner bowtie2)``). Python packages are
printed as ``[OK]``/``[MISSING]`` with their name only. With ``-q`` the tool lines also show only
the status and name, and the section headings and notes (logged at INFO) are hidden.

Missing optional items are summarised in a note and do not fail the check. The command exits with
an error (status 2) only when a tool or package from the default-pipeline group is missing.

If ``bowtie2`` is the conda build -- detected by the absence of the AVX2 binary
``bowtie2-align-s-v256`` next to ``bowtie2`` -- a note suggests rebuilding it from source for
AVX2 speed: ``pixi run install-bowtie2`` in a pixi environment, or
``bash install_scripts/pixi_install_bowtie2.sh`` in a conda environment.

Invocation
==========

.. code-block:: text

    AAFTF dependency [-q] [-v]

Example
=======

.. code-block:: bash

    # full report, with tool locations / where missing tools are used
    AAFTF dependency

    # only the OK/MISSING lines; the exit status tells whether the default pipeline can run
    AAFTF dependency -q && echo "ready"
