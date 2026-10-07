========
database
========

Lists and fetches the reference databases AAFTF's other subcommands rely on, into the database
folder (``$AAFTF_DB``, or ``~/.cache/aaftf`` when it is unset). This is the only command that
downloads them: ``filter``/``vecscreen``/``taxpurge``/``fcs_screen`` read the stored copies and,
if one is missing, stop with the ``AAFTF database ...`` command that fetches it. (The one
exception is ``filter``'s ``-a/--screen_accessions`` and ``-u/--screen_urls``, which fetch those
extra sequences themselves.)

Where databases are stored
==========================

* ``$AAFTF_DB`` unset: databases go to ``~/.cache/aaftf`` (or ``$XDG_CACHE_HOME/aaftf``). Some are
  large -- each sourmash index and the FCS container image are several GB -- so ``database``
  warns before filling your home directory, which on clusters often has a small quota.
* ``$AAFTF_DB`` may list several folders separated by ``:``, like ``$PATH``. Each database is read
  from the first folder that has it, and missing ones are downloaded into the first folder you can
  write to. This lets a lab share one read-only copy while each user keeps a personal folder for
  anything missing::

      export AAFTF_DB=/shared/lab/aaftf_db:/scratch/$USER/aaftf_db

* If none of the listed folders is writable, downloads fall back to ``~/.cache/aaftf`` (with the
  same warning).


Databases
=========

Run ``AAFTF database`` with no arguments to list every database: its short name, file, size,
which subcommands use it, and the folder it is stored in (or ``not downloaded``, with the
remote size). Like ``AAFTF dependency``, the list is split in two: the databases ``AAFTF pipeline_short``
needs with its default settings (``phix``, ``univec``, ``euks``, ``proks``, ``mitodb``),
then those only needed by non-default options or optional steps (``sm_gbk``, ``sm_gtdbrep``, ``sm_gtdb``,
``fcs_script``, ``fcs_image``); any other files found in the database folders are listed last.

.. list-table::
   :header-rows: 1
   :widths: 14 40 46

   * - Name
     - File
     - Used by
   * - ``phix``
     - ``GCF_000819615.1_ViralProj14015_genomic.fna.gz`` (PhiX genome)
     - ``filter``
   * - ``univec``
     - ``UniVec``
     - ``filter``, ``vecscreen``
   * - ``euks``
     - ``contam_in_euks.fa.gz``
     - ``vecscreen``
   * - ``proks``
     - ``contam_in_prok.fa``
     - ``vecscreen``
   * - ``mitodb``
     - ``mitochondrion.1.1.genomic.fna.gz`` (RefSeq mitochondria)
     - ``vecscreen``
   * - ``sm_gbk``
     - ``genbank-k31.lca.json.gz`` (sourmash GenBank)
     - ``taxpurge --sourdb_type gbk``
   * - ``sm_gtdbrep``
     - ``gtdb-rs220-reps.k31.lca.json.gz`` (sourmash GTDB representatives)
     - ``taxpurge --sourdb_type gtdbrep``
   * - ``sm_gtdb``
     - ``gtdb-rs220-k31.lca.json.gz`` (sourmash GTDB)
     - ``taxpurge --sourdb_type gtdb``
   * - ``fcs_script``
     - ``run_fcsadaptor.sh``
     - ``fcs_screen``
   * - ``fcs_image``
     - ``fcs-adaptor.0.5.5.sif``
     - ``fcs_screen`` (singularity)

The sourmash databases and the FCS image are several GB each; download only the sourmash index
matching how you'll call ``taxpurge --sourdb_type``.

Downloads are safe to interrupt: each file is written to a ``.tmp`` path and renamed on success,
and files already present in any database folder are skipped on re-run unless ``--force`` is
given.

Invocation
==========

.. code-block:: text

    AAFTF database [DATABASE ...] [--force] [-q] [-v]

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Description
   * - ``DATABASE``
     - Databases to download, by short name (e.g. ``univec``) or file name (e.g. ``UniVec``), case
       insensitive; or ``required`` (what the default pipeline needs), ``optional`` (the rest) or
       ``all``. With none, list the databases instead.
   * - ``--force``
     - Re-download even if already present (into the first writable folder).

Example
=======

.. code-block:: bash

    # Pick a folder with plenty of space (optional; defaults to ~/.cache/aaftf)
    export AAFTF_DB=/path/with/space/aaftf_db

    # See what is available, what is already downloaded, and where
    AAFTF database

    # Everything the default pipeline needs (phix univec euks proks mitodb)
    AAFTF database required

    # ... plus the sourmash index taxpurge uses (pipelines run with -p/--phylum)
    AAFTF database required sm_gbk

    # ... plus the FCS-adaptor files for fcs_screen
    AAFTF database required fcs_script fcs_image

    # Names and file names can be mixed
    AAFTF database UniVec contam_in_prok.fa

    # Everything (all sourmash indices + FCS-adaptor); a lot of disk
    AAFTF database all

Container usage
================

The Singularity definition (``AAFTF.def``) defaults to ``skip_db_download=1``, so the image
contains no databases. Build with ``--build-arg skip_db_download=0`` to run ``AAFTF database phix
univec euks proks mitodb sm_gbk fcs_script fcs_image`` during ``%post`` and bake them into the
image at ``/opt/aaftf_db``. Otherwise (and for the Docker image, which never includes them),
bind-mount a database folder at ``/opt/aaftf_db`` (the container's ``$AAFTF_DB``) and run
``AAFTF database NAME ...`` into it before using ``filter``/``vecscreen``/``taxpurge``.
