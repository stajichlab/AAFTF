============
fcs_gx_purge
============

Purges contaminant contigs from one or more assemblies using NCBI's **FCS-GX** tool -- a genome
cross-contamination screen that classifies every contig/region against a comprehensive reference
database of taxonomy-labeled genomes, rather than the k-mer/sketch-based approach used by
:doc:`taxpurge`.

FCS-GX needs its whole database (~470 GB for the full ``gxdb``) in RAM, so ``fcs_gx_purge`` is a
standalone utility (listed under "Utilities" in ``AAFTF --help``), **not** a pipeline step: the
``pipeline_*`` subcommands use :doc:`taxpurge` (with ``-p``) as their taxonomic screen. Run
``fcs_gx_purge`` separately, on a large-memory node, on finished assemblies -- ideally several at
once, since the database is set up only once per run.

Algorithm
=========

1. Checks that ``run_gx.py`` is on ``PATH``, that every input exists and that the database's
   ``.gxi`` exists (``FileNotFoundError``, exit code 2, otherwise) -- before any time is spent
   copying the database.
2. Unless the database is on tmpfs (e.g. ``/dev/shm``, where it is already in RAM), checks that
   the job has enough RAM to keep the whole database cached: its size plus 5%. Available RAM
   respects SLURM / cgroup memory limits. If it is too small, the run stops before anything is
   copied, with a ``MemoryError`` (exit code 1) that suggests a memory request, e.g.
   ``sbatch --mem=...G``.
3. Decides where ``run_gx.py`` reads the database from:

   * on tmpfs (``/dev/shm``): used in place;
   * otherwise a quick **speed test** reads random 4 KiB blocks from the largest database file
     (the ``.gxs``) -- FCS-GX's own access pattern -- for at most 256 reads or 3 seconds. A mean
     read time of at most 0.2 ms (local NVMe/SSD, or data already in the page cache) means the
     database is used in place;
   * slower storage (network filesystems such as Lustre/GPFS/NFS, spinning disks) is copied once
     -- every file sharing the database prefix -- into a new per-job folder
     ``aaftf-gxdb-<SLURM_JOB_ID>-*`` in ``--copy_to``: by default ``$TMPDIR`` (under SLURM
     usually the job's own node-local folder, which SLURM removes at the end of the job), or
     ``/scratch`` when ``$TMPDIR`` is not set. Its free space is checked first (``OSError``
     otherwise). The copy is skipped, with a warning, when ``--copy_to`` is on the same
     filesystem as the database (it could not be faster), and made with a warning when
     ``--copy_to`` is on a network filesystem.

   The log reports the measured read time and filesystem type.
4. For each assembly, runs ``run_gx.py --fasta INPUT --tax-id TAXID --gx-db DB --out-dir DIR``
   with ``GX_NUM_CORES`` set to ``-c/--cpus``, and parses its report (tab-separated: ``seq_id start_pos end_pos seq_len
   action div agg_cont_cov top_tax_name``). Only rows whose action is ``EXCLUDE`` mark a sequence
   for removal; ``TRIM``, ``FIX``, ``REVIEW`` and other actions are logged as warnings and the
   sequence is **kept** (partial-contig edits are not applied -- review them manually).
5. Removes the copy at the end, also when a run fails. A job killed outright (out of memory,
   time limit) cannot clean up; the default ``$TMPDIR`` is removed by SLURM on most clusters.
   Otherwise remove leftover ``aaftf-gxdb-*`` folders yourself.

Each assembly's ``run_gx.py`` stderr goes to ``fcs_gx.log`` in its own subfolder of the work
directory, which also holds ``fcs_gx_purge.log`` and is kept when given with ``-w``, with ``-v``,
or when the run fails. If ``run_gx.py`` exits non-zero or writes no report, the command fails with
a ``RuntimeError`` (exit code 1).

Outputs
=======

For each input, with ``<prefix>`` its file name without ``.fasta``/``.fa``/``.fna``/``.fas``/``.fsa``
(and an optional ``.gz``), written to ``-o/--outdir`` or, by default, next to the input:

* ``<prefix>.fcsgx_cleaned.fasta`` -- the assembly without the ``EXCLUDE`` contigs
* ``<prefix>.fcs_gx_report.txt`` -- the FCS-GX report, for later review

Two inputs with the same ``<prefix>`` going to the same folder are an error (``ValueError``).

Requires the FCS-GX database
============================

FCS-GX needs a very large reference database set up separately per NCBI's instructions
(https://github.com/ncbi/fcs/wiki/FCS-GX) -- AAFTF does not download or manage this database for
you. ``-d/--db`` (required) is the database path prefix, e.g. ``/path/gxdb/all`` for
``/path/gxdb/all.gxi``, ``all.gxs``, ``all.taxa.tsv``, ...; the ``.gxi`` file itself or the folder
holding it (if it has only one ``.gxi``) also work. The database is all those files together, not
only the ``.gxi``: they are measured, tested and copied as one. Plan for a job memory request of about the database size plus 5% plus a
little headroom (roughly 500 GB for the full ``gxdb``).

Cutoffs / defaults
==================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-d/--db``
     - (required)
     - FCS-GX database path prefix
   * - ``-t/--taxid``
     - (required)
     - NCBI Taxonomy ID for the query organism or a parent taxon (e.g. 4890 = Ascomycota); FCS-GX
       uses this to decide what counts as "self" vs. contaminant, so a wrong ID can drop real contigs
   * - ``-c/--cpus``
     - CPUs allocated to the job
     - Threads for ``run_gx.py`` (``GX_NUM_CORES``)
   * - ``-o/--outdir``
     - next to each input
     - Directory for the cleaned assemblies and FCS-GX reports
   * - ``--copy_to``
     - ``$TMPDIR``, else ``/scratch``
     - Fast node-local folder for the per-job database copy
   * - Speed test
     - 0.2 ms
     - Largest mean random 4 KiB read time at which the database is used in place
   * - Memory check
     - database size + 5%
     - RAM available to the job (not checked for a database on tmpfs)
   * - Contig removal criterion
     - ``EXCLUDE`` action
     - FCS-GX decides the ``action`` internally; AAFTF drops only ``EXCLUDE`` sequences and
       warns about (but keeps) everything else

Invocation
==========

.. code-block:: text

    AAFTF fcs_gx_purge -i FASTA [FASTA ...] -d DB -t TAXID [-o OUTDIR] [-w DIR] [-c INT]
                       [--copy_to DIR] [-q] [-v]

``-i/--input`` (one or more assemblies), ``-d/--db`` and ``-t/--taxid`` are required.

Example
=======

A SLURM job screening several finished assemblies; a database on slow shared storage is copied
once to the job's own ``$TMPDIR``:

.. code-block:: bash

    #!/bin/bash -l
    #SBATCH --mem=500G -c 16 --time=12:00:00 -J fcsgx

    AAFTF fcs_gx_purge -d /shared/db/gxdb/all -t 4890 -o genomes/fcsgx \
        -i genomes/STRAINX.final.fasta genomes/STRAINY.final.fasta genomes/STRAINZ.final.fasta

    # writes genomes/fcsgx/STRAINX.fcsgx_cleaned.fasta and genomes/fcsgx/STRAINX.fcs_gx_report.txt,
    # and the same for STRAINY and STRAINZ

Several rounds with different taxonomy IDs
==========================================

``-t/--taxid`` applies to every assembly of a run, so assemblies of different organisms need
separate runs. A database copy that ``fcs_gx_purge`` makes itself is removed when each run ends,
so every round would copy the ~470 GB database again. Instead, copy it once yourself to fast
storage -- ``/dev/shm`` or a node-local ``$TMPDIR`` -- and pass that copy with ``-d``:
``fcs_gx_purge`` uses a database it did not copy in place and never removes it. Remove the copy
yourself when you are done:

.. code-block:: bash

    #!/bin/bash -l
    #SBATCH --mem=500G -c 16 --time=24:00:00 -J fcsgx

    GXDB=$TMPDIR/gxdb            # or /dev/shm/gxdb (needs --mem for the database too)
    mkdir -p "$GXDB"
    cp /shared/db/gxdb/all.* "$GXDB"/
    trap 'rm -rf "$GXDB"' EXIT    # remove the copy when the job script ends

    AAFTF fcs_gx_purge -d "$GXDB/all" -t 4890 -o genomes/fcsgx -i genomes/ASCO*.final.fasta
    AAFTF fcs_gx_purge -d "$GXDB/all" -t 5204 -o genomes/fcsgx -i genomes/BASIDIO*.final.fasta

A copy in ``/dev/shm`` stays in RAM if the job is killed before ``trap`` runs; a copy in the job's
``$TMPDIR`` is removed by SLURM on most clusters.
