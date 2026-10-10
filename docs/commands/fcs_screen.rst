==========
fcs_screen
==========

Runs NCBI's official **FCS-adaptor** tool (the same tool NCBI uses during genome submission QC)
to detect and trim residual sequencing-adaptor/vector contamination from assembled contigs. This
is an alternative (NCBI-authoritative) take on part of what :doc:`vecscreen` does, delegated to
NCBI's containerized tool rather than reimplemented in AAFTF.

Algorithm
=========

AAFTF is a thin wrapper. The ``run_fcsadaptor.sh`` launcher is taken from ``--fcs_script``, else
from ``$PATH``, else from the ``fcs_script`` database installed by ``AAFTF database``; the
Singularity image likewise comes from ``--image`` or the ``fcs_image`` database (a missing
database stops the run with exit code 2). It then invokes::

    run_fcsadaptor.sh --fasta-input INFILE --output-dir WORKDIR \
        {--euk|--prok} --container-engine {singularity|docker} --image IMAGE

``--euk`` (eukaryotic screening) is used unless ``--prok`` is passed. The container itself
downloads nothing further at run time -- it ships a self-contained adaptor database. AAFTF then
copies FCS-adaptor's ``cleaned_sequences/<input basename>`` to ``-o/--outfile``, and renames its
``fcs_adaptor_report.txt`` to ``{outfile}.fcs_adaptor_report.txt`` (printed to stdout as well).

Errors: a missing container engine raises ``FileNotFoundError`` (exit code 2), an unknown
``--container_engine`` a ``ValueError``, and an FCS-adaptor failure (non-zero exit, or no cleaned
FASTA written) a ``RuntimeError`` (exit code 1; rerun with ``-v`` to see the tool output). The
work directory holds ``fcs_screen.log`` and is kept when given with ``-w``, with ``-v``, or when
the run fails.

Requires a container engine
=============================

Unlike every other AAFTF subcommand, ``fcs_screen`` requires **an additional container runtime**
on the host (or inside the AAFTF container, when using Docker-in-Docker or Singularity nested
execution) -- NCBI ships FCS-adaptor only as a container image:

* ``--container_engine singularity`` (default): requires ``singularity`` or ``apptainer`` on
  ``$PATH``. The ``.sif`` image is the ``fcs_image`` database from ``AAFTF database`` (or
  ``--image PATH``).
* ``--container_engine docker``: requires ``docker`` on ``$PATH``. AAFTF passes a registry
  reference (``ncbi/fcs-adaptor:{VERSION}``, unless ``--image`` is given) for Docker to
  pull/cache itself.

Pinned tool version: FCS-adaptor **0.5.5** (hard-coded in ``aaftf/resources.py``).

Invocation
==========

.. code-block:: text

    AAFTF fcs_screen -i FASTA -o FASTA [-w DIR] [--image IMAGE]
                     [--container_engine {singularity,docker}] [--prok]
                     [--fcs_script FCS_SCRIPT] [-q] [-v]

``-i/--input`` (alias ``--infile``) and ``-o/--outfile`` are required. Eukaryotic mode is the
default; ``--prok`` switches to prokaryotic screening.

Example
=======

.. code-block:: bash

    # inside a SLURM/HPC node with singularity available
    AAFTF fcs_screen -i genomes/STRAINX.spades.fasta -o genomes/STRAINX.fcs_adaptor.fasta

    # docker instead
    AAFTF fcs_screen --container_engine docker \
        -i genomes/STRAINX.spades.fasta -o genomes/STRAINX.fcs_adaptor.fasta

.. note::
   ``fcs_screen`` is not currently wired into ``AAFTF pipeline_short`` -- run it as an extra manual step
   (typically in place of, or alongside, :doc:`vecscreen`) if you need NCBI-authoritative adaptor
   screening prior to genome submission. See also :doc:`fcs_gx_purge` for NCBI's genomic
   cross-contamination screen, and :doc:`fix_tbl` for adjusting NCBI ``.tbl`` annotation
   coordinates after FCS trims a sequence.
