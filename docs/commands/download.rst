========
download
========

**Aliases:** ``configure``, ``install``, ``download_db``, ``setup``

Fetches and caches every reference database AAFTF's other subcommands rely on into a persistent
``$AAFTF_DB`` directory, so ``filter``/``vecscreen``/``sourpurge``/``fcs_screen`` don't each
re-download (large) files the first time they run.

What it downloads
==================

Unless skipped with the flags below, ``download`` fetches:

* **Core contamination databases** (``--skip-core`` to disable): PhiX genome, UniVec, NCBI
  ``contam_in_euks``/``contam_in_prok``, RefSeq mitochondrion sequences -- used by ``filter`` and
  ``vecscreen``.
* **sourmash LCA taxonomy databases** (``--skip-sourmash`` to disable): one of GenBank k=31
  (``gbk``), GTDB (``gtdb``), the smaller GTDB representative-genomes set (``gtdbrep``), or
  ``all`` -- used by ``sourpurge``. These are large (GenBank is several GB); pick the one that
  matches how you'll call ``sourpurge --sourdb_type``.
* **NCBI FCS-adaptor resources** (``--skip-fcs`` to disable): the ``run_fcsadaptor.sh`` wrapper
  script plus a cached Singularity ``.sif`` image -- used by ``fcs_screen``.

Downloads are resumable/idempotent: each file is written to a ``.tmp`` path and atomically
renamed on success, and files that already exist are skipped on re-run unless ``--force`` is
given.

Invocation
==========

.. code-block:: text

    AAFTF download [--AAFTF_DB DIR] [--force] [--skip-core] [--skip-sourmash]
                   [--sourdb-type {gbk,gtdb,gtdbrep,all}] [--skip-fcs] [-v] [--pipe]

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Description
   * - ``--AAFTF_DB DIR``
     - Target directory. Defaults to the ``$AAFTF_DB`` environment variable.
   * - ``--force``
     - Re-download files even if already present.
   * - ``--skip-core``
     - Skip UniVec/PhiX/Euk/Prok/Mito contamination databases.
   * - ``--skip-sourmash``
     - Skip sourmash taxonomy databases.
   * - ``--sourdb-type``
     - ``gbk`` (default here is ``all``; ``sourpurge`` itself defaults to ``gbk``), ``gtdb``,
       ``gtdbrep``, or ``all``.
   * - ``--skip-fcs``
     - Skip NCBI FCS-adaptor script + container image.

Example
=======

.. code-block:: bash

    export AAFTF_DB=~/lib/AAFTF_DB
    mkdir -p "$AAFTF_DB"

    # Everything (core + all sourmash indices + FCS-adaptor); can take a while / a lot of disk
    AAFTF download --AAFTF_DB "$AAFTF_DB"

    # Just the databases needed for filter/vecscreen/sourpurge with the GenBank sourmash index
    AAFTF download --AAFTF_DB "$AAFTF_DB" --sourdb-type gbk --skip-fcs

Container usage
================

The Singularity build (``AAFTF.def``) calls this exact command during ``%post`` (unless built with
``--build-arg skip_db_download=1``) to bake the GenBank sourmash database into the image at
``/opt/aaftf_db``. If you build with ``skip_db_download=1`` (the CI default, to keep the image
small), run ``AAFTF download`` yourself against a bind-mounted directory before using
``sourpurge``/``vecscreen``/``filter`` from that image.
