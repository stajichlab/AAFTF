========
database
========

Lists and fetches the reference databases AAFTF's other subcommands rely on, into the database
folder (``$AAFTF_DB``, or ``~/.cache/aaftf`` when it is unset). This is the only command that
downloads them: ``filter``/``vecscreen``/``taxpurge``/``fcs_screen`` read the stored copies and,
if one is missing, stop with the ``AAFTF database ...`` command that fetches it. (The one
exception is ``filter``'s ``-a/--screen_accessions`` and ``-u/--screen_urls``, which fetch those
extra sequences themselves.)

The contamination database (``contamdb``)
==========================================

``filter`` (reads, before assembly) and ``vecscreen`` (contigs, after assembly) screen against the
same database. ``AAFTF database contamdb`` (also part of ``required``) builds it, rather than
downloading it, as ``contamdb.fa`` plus a BLAST database next to it (``contamdb.nhr``, ``.nin``,
``.nsq``, ...), from:

1. NCBI's adaptor set for screening eukaryote assemblies, ``adaptors_for_screening_euks.fa``
   (downloaded as ``euk_adaptors`` if missing): about 600 Illumina, 454 and other sequencing
   adaptors, primers and barcodes, 20-70 bp each.
2. The PhiX genome (downloaded as ``phix`` if missing).
3. The BBTools ``lambda.fa.gz`` (lambda phage spike-in) and ``adapters.fa``, from the
   ``resources/`` folder next to ``bbduk.sh``.

Both steps only require ``contamdb.fa``, so you can add your own sequences to it. ``vecscreen``
checks the BLAST database before it starts and rebuilds it with ``makeblastdb`` when a file is
missing, empty or older than ``contamdb.fa``. ``AAFTF database contamdb`` keeps an existing
``contamdb.fa`` (it only updates the BLAST database); ``--force`` rebuilds it from its sources,
discarding any sequences you added. Duplicate and contained sequences are removed with BBTools
``dedupe.sh``.

``contamdb`` holds no cloning vectors. UniVec's vectors carry yeast and other fungal DNA
(ADH1/GAL4 parts, LEU2, TRP1, HIS3, 2-micron, CEN/ARS, telomere seeds in YAC/TAR/two-hybrid
vectors) under names that don't say so; screening against them split genuine fungal contigs (the
S. cerevisiae S288C genome came out of ``vecscreen`` as 150 pieces instead of 17 chromosomes) and
removed genuine fungal reads in ``filter``. To screen for a particular vector anyway, add its
sequence to ``contamdb.fa``.

Once ``contamdb`` is built, the downloaded NCBI adaptor and PhiX files are deleted (only from the
folder AAFTF downloads into; a copy in another ``$AAFTF_DB`` folder may be shared and is kept).
``AAFTF database`` then lists them as ``fulfilled by contamdb`` and skips them when asked for again
(``required``, ``all``, or by name), just as a sourmash sketch file is skipped once its RocksDB
index is built (``fulfilled by rocksdb``, also with ``--no-build``); ``--force`` downloads them.

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
needs with its default settings (``phix`` and ``euk_adaptors``, which ``contamdb`` is built from --
shown as ``fulfilled by contamdb`` once it is built, even if their files were removed -- then
``contamdb``), then the rest (the sourmash databases and their lineage CSVs,
``fcs_script``, ``fcs_image``); anything else directly in the database folders is listed last (a subfolder as one
entry, without its contents or size).

.. list-table::
   :header-rows: 1
   :widths: 14 40 46

   * - Name
     - File
     - Used by
   * - ``phix``
     - ``GCF_000819615.1_ViralProj14015_genomic.fna.gz`` (PhiX genome)
     - ``contamdb`` (source)
   * - ``euk_adaptors``
     - ``adaptors_for_screening_euks.fa`` (NCBI adaptors for eukaryote assemblies)
     - ``contamdb`` (source)
   * - ``contamdb``
     - ``contamdb.fa`` and its BLAST database ``contamdb``, built locally (see above)
     - ``filter``, ``vecscreen``
   * - ``gbk-fungi`` (+ ``gbk-fungi-lineages``)
     - ``genbank-20250408-fungi-k31.zip`` (5.3 GB; GenBank 2025.04 fungi, k=31, scaled=1000) and
       ``genbank-20250408-lineages.fungi.csv``
     - ``taxpurge --sourdb gbk-fungi``
   * - ``gbk-fungi-lite`` (+ ``gbk-euks-lineages``)
     - ``genbank-euks-2025.01-fungi.k31.sig.zip`` (170 MB; GenBank 2025.01 fungi, scaled=10000) and
       ``genbank-euks-2025.01-eukaryotes.lineages.csv``
     - ``taxpurge --sourdb gbk-fungi-lite`` (default)
   * - ``gtdb-reps`` (+ ``gtdb-reps-lineages``)
     - ``gtdb-rs226-reps.k31.sig.zip`` (3.9 GB; GTDB rs226 representatives, bacteria and archaea,
       scaled=1000) and ``gtdb-rs226-reps.lineages.csv``
     - ``taxpurge --sourdb gtdb-reps``
   * - ``gtdb-reps-lite`` (+ ``gtdb-reps-lineages``)
     - ``gtdb-rs226-reps.k31-sc10k.sig.zip`` (447 MB; the same genomes at scaled=10000)
     - ``taxpurge --sourdb gtdb-reps-lite`` (default)
   * - ``human-hg38``
     - ``hg38-entire.sig.zip`` (55 MB; human hg38, scaled=1000)
     - ``taxpurge --sourdb human-hg38`` (default)
   * - ``fcs_script``
     - ``run_fcsadaptor.sh``
     - ``fcs_screen``
   * - ``fcs_image``
     - ``fcs-adaptor.0.5.5.sif``
     - ``fcs_screen`` (singularity)

The sourmash databases and the FCS image are several GB each. ``taxpurge`` needs both files
(sketches and lineages) of each database it searches; naming a sourmash database downloads its
lineage CSV with it, and ``taxpurge`` downloads the default set (``gbk-fungi-lite``,
``gtdb-reps-lite``, ``human-hg38``; ~0.7 GB). Sourmash databases are installed in their final
form, a RocksDB index (searched from disk by ``taxpurge``, with much less memory): ``gtdb-reps`` and
``gtdb-reps-lite`` as the published index (its archive is downloaded and unpacked in a temporary
folder, so an interrupted run leaves nothing half-done), the others by downloading the sketch file,
building the index from it (``sourmash scripts index``) and removing the sketch file, which is then
listed as "fulfilled by rocksdb", with the index on the row below it. ``--no-build`` downloads
those sketch files without building their indexes (``taxpurge`` then searches the sketch files).
See :doc:`taxpurge`.

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
     - Databases to download, by short name (e.g. ``contamdb``) or file name (e.g. ``contamdb.fa``), case
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

    # Everything the default pipeline needs (phix euk_adaptors contamdb)
    AAFTF database required

    # ... plus the sourmash databases taxpurge searches (pipelines run with -p/--phylum)
    AAFTF database required taxpurge

    # ... plus the FCS-adaptor files for fcs_screen
    AAFTF database required fcs_script fcs_image

    # Names and file names can be mixed
    AAFTF database contamdb.fa gtdb-reps-lite

    # Everything (all sourmash indices + FCS-adaptor); a lot of disk
    AAFTF database all

Container usage
================

The Singularity definition (``AAFTF.def``) defaults to ``skip_db_download=1``, so the image
contains no databases. Build with ``--build-arg skip_db_download=0`` to run ``AAFTF database required
taxpurge fcs_script fcs_image`` during ``%post`` and bake them into the
image at ``/opt/aaftf_db``. Otherwise (and for the Docker image, which never includes them),
bind-mount a database folder at ``/opt/aaftf_db`` (the container's ``$AAFTF_DB``) and run
``AAFTF database NAME ...`` into it before using ``filter``/``vecscreen``/``taxpurge``.
