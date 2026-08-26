==========
fcs_screen
==========

**Aliases:** ``ncbi_fcs``, ``ncbi_fcs-screen``

Runs NCBI's official **FCS-adaptor** tool (the same tool NCBI uses during genome submission QC)
to detect and trim residual sequencing-adaptor/vector contamination from assembled contigs. This
is an alternative (NCBI-authoritative) take on part of what :doc:`vecscreen` does, delegated to
NCBI's containerized tool rather than reimplemented in AAFTF.

Algorithm
=========

AAFTF is a thin wrapper: it locates or downloads the ``run_fcsadaptor.sh`` launcher script (from
``$AAFTF_DB`` or NCBI's GitHub release, cached under ``$AAFTF_DB``) and a container image, then
invokes::

    run_fcsadaptor.sh --fasta-input INFILE --output-dir WORKDIR \
        {--euk|--prok} --container-engine {singularity|docker} --image IMAGE

``--euk`` (eukaryotic screening) is used unless ``--prok`` is passed. The container itself
downloads nothing further at run time -- it ships a self-contained adaptor database. AAFTF then
copies FCS-adaptor's ``cleaned_sequences/<input basename>`` to ``-o/--outfile``, and renames its
``fcs_adaptor_report.txt`` to ``{outfile}.fcs_adaptor_report.txt`` (printed to stdout as well).

Requires a container engine
=============================

Unlike every other AAFTF subcommand, ``fcs_screen`` requires **an additional container runtime**
on the host (or inside the AAFTF container, when using Docker-in-Docker or Singularity nested
execution) -- NCBI ships FCS-adaptor only as a container image:

* ``--container_engine singularity`` (default): requires ``singularity`` or ``apptainer`` on
  ``$PATH``. The ``.sif`` image is downloaded once to
  ``$AAFTF_DB/fcs-adaptor.{VERSION}.sif`` and reused thereafter.
* ``--container_engine docker``: requires ``docker`` on ``$PATH``. AAFTF passes a registry
  reference (``ncbi/fcs-adaptor:{VERSION}``) for Docker to pull/cache itself.

Pinned tool version: FCS-adaptor **0.5.5** (hard-coded in ``AAFTF/resources.py``).

Invocation
==========

.. code-block:: text

    AAFTF fcs_screen -i INFILE -o OUTFILE
                     [--container_engine {singularity,docker}] [--image IMAGE]
                     [--euk | --prok] [--fcs_script PATH]
                     [--AAFTF_DB DIR] [-w WORKDIR] [-v] [--pipe]

``-i/--input`` and ``-o/--outfile`` are required.

Example
=======

.. code-block:: bash

    # inside a SLURM/HPC node with singularity available
    AAFTF fcs_screen -i genomes/STRAINX.spades.fasta -o genomes/STRAINX.fcs_adaptor.fasta \
        --AAFTF_DB "$AAFTF_DB"

    # docker instead
    AAFTF fcs_screen --container_engine docker \
        -i genomes/STRAINX.spades.fasta -o genomes/STRAINX.fcs_adaptor.fasta

.. note::
   ``fcs_screen`` is not currently wired into ``AAFTF pipeline`` -- run it as an extra manual step
   (typically in place of, or alongside, :doc:`vecscreen`) if you need NCBI-authoritative adaptor
   screening prior to genome submission. See also :doc:`fcs_gx_purge` for NCBI's genomic
   cross-contamination screen, and :doc:`fix_tbl` for adjusting NCBI ``.tbl`` annotation
   coordinates after FCS trims a sequence.
