============
Installation
============

AAFTF can be run three ways: a **conda/pip environment** you assemble yourself, a **pixi**-managed
environment (the same mechanism used to build the project's own Docker/Singularity images), or
directly from a pre-built **Docker/Singularity/Apptainer container**. All external bioinformatics
tools (SPAdes, BBTools, BLAST+, sourmash, bwa, minimap2, pypolca, polypolish, NextPolish2, mosdepth, ...)
must be reachable on ``$PATH`` -- the three methods below differ only in *how* that PATH gets populated.

Option 1: conda environment
============================

.. code-block:: bash

    # only what the default pipeline needs (same as the pixi "default" environment);
    # installs AAFTF from PyPI
    conda env create -f environment.yml
    conda activate aaftf

    # or everything AAFTF can use (same as the pixi "complete" environment), with AAFTF
    # installed from this checkout in editable mode
    conda env create -f environment.dev.yml
    conda activate aaftf-dev

.. warning::
   AAFTF requires samtools >= 1.13 (the environment files pin >= 1.24). If you assemble your own
   environment, watch for bioconda packages that bundle an old, conflicting samtools binary of
   their own (older MaSuRCA builds were a historical example of this) -- AAFTF's own dependency
   lists (``pyproject.toml``, ``environment.yml``) no longer pull in any such package.

Option 2: pixi (recommended for reproducible/locked environments)
====================================================================

The repository ships a ``pyproject.toml`` (``[tool.pixi.*]`` tables) / ``pixi.lock`` pinning every dependency (this is what the
Docker and Singularity images are built from). From a checkout of the repository:

.. code-block:: bash

    curl -fsSL https://pixi.sh/install.sh | bash   # install pixi, if needed
    cd AAFTF
    pixi install -e complete               # every tool AAFTF can use
    pixi shell -e complete                 # activate the environment
    AAFTF --version

Three environments are defined, each installing AAFTF from the checkout in editable mode:

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Environment
     - Features
     - Contents
   * - ``default``
     - default
     - Only what ``AAFTF pipeline_short``, ``pipeline_long`` and ``pipeline_hybrid`` need with their
       default settings (BBTools, SPAdes, BLAST+, bwa, samtools, minimap2, Flye, Racon,
       Polypolish, pypolca, freebayes). ``pixi install`` alone creates this one.
   * - ``complete``
     - default, optional
     - Plus every tool used by optional steps (``polish``, ``depth``, ``mito``, ``fcs_screen``,
       ``fcs_gx_purge``) and non-default options (fastp, Trimmomatic, MEGAHIT, Unicycler, ...).
       The Docker and Singularity images use this environment.
   * - ``dev``
     - default, optional, dev
     - Plus pytest, pytest-cov and pre-commit, for working on AAFTF itself.

``AAFTF dependency`` reports what an environment is missing.

The ``complete`` and ``dev`` environments include bowtie2 (used by ``filter --aligner bowtie2``) from
conda. That build cannot use AVX2 instructions (it falls back from its x86-64-v3 version at
runtime), so it runs slower than it could. If you use bowtie2 and want AVX2 support, rebuild it from
source into the environment once (this replaces the conda binaries; the Docker and Singularity
images already do this):

.. code-block:: bash

    pixi run -e complete install-bowtie2   # or: pixi run -e dev install-bowtie2

Option 3: Docker / Singularity / Apptainer
============================================

Prebuilt containers bundle every dependency, including ``pypolca`` (POLCA-style polishing that
works with modern samtools) and a source-built ``bowtie2``. Only ``fcs_screen`` additionally requires a
*separate* container engine (singularity/apptainer or docker) at runtime, since NCBI FCS-adaptor
itself ships as a container image invoked from inside AAFTF -- see :doc:`commands/fcs_screen`.

**Docker**

.. code-block:: bash

    # Build the pixi "complete" environment (AAFTF_VERSION must be supplied explicitly, as .git is
    # excluded from the build context; it is converted to a PEP 440 version for the package).
    # Add --build-arg PIXI_ENV=default for a smaller image with only the default-pipeline tools.
    docker build --build-arg AAFTF_VERSION=$(git describe --tags --always) -t aaftf:latest .

    # Run
    docker run --rm aaftf:latest AAFTF --help

    # Run with data + a persistent reference-database directory mounted in
    docker run --rm \
        -v /path/to/data:/data \
        -v /path/to/aaftf_db:/opt/aaftf_db \
        aaftf:latest AAFTF trim --read1 /data/R1.fq.gz --read2 /data/R2.fq.gz

A second, simpler conda-based image is available via ``Dockerfile.conda`` for environments where
the pixi-based build is undesirable; build/run commands are the same.

**Singularity / Apptainer**

.. code-block:: bash

    # Build (defaults: git_ref=main, skip_db_download=1 -- no databases in the image; bind-mount
    # or download AAFTF_DB separately at run time, or add --build-arg skip_db_download=0 to bake them in)
    singularity build AAFTF.sif AAFTF.def

    # Build a specific tagged release
    singularity build --build-arg git_ref=v0.7.0 --build-arg skip_db_download=1 AAFTF.sif AAFTF.def

    # Run
    singularity exec AAFTF.sif AAFTF --help
    singularity exec --bind /path/to/data:/data AAFTF.sif \
        AAFTF trim --read1 /data/R1.fq.gz --read2 /data/R2.fq.gz

    # Or use the built-in runscript
    singularity run --bind /path/to/data:/data AAFTF.sif trim --read1 /data/R1.fq.gz --read2 /data/R2.fq.gz

The database directory defaults to ``/opt/aaftf_db`` inside the container (override with
``AAFTF_DB``). Both the Docker entrypoint and the Singularity ``%environment``/``/etc/profile.d``
hook ensure the pixi-managed tool PATH survives both interactive (``singularity exec``) and login
shells (``bash -l``, as used by SLURM/Nextflow task scripts).

External tools
==============

**Needed by the pipelines** (``pipeline_short``, ``pipeline_long`` and ``pipeline_hybrid`` with their
default settings):

* `BBTools <https://bbmap.org/>`_ (``bbduk.sh``, ``dedupe.sh``, ``reformat.sh``; needs Java) -- read
  trimming and filtering, building the ``contamdb`` database
* `SPAdes <https://github.com/ablab/spades>`_ -- Illumina assembly
* `NCBI BLAST+ <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ -- adaptor/vector
  screening (``vecscreen``)
* `bwa <https://github.com/lh3/bwa>`_ and `samtools <https://github.com/samtools/samtools>`_ >= 1.13
  -- read mapping for ``covpurge``
* `minimap2 <https://github.com/lh3/minimap2>`_ -- duplicate contig detection (``rmdup``) and read
  mapping for ``depth``
* `mosdepth <https://github.com/brentp/mosdepth>`_ -- read depth for ``depth``, which the pipelines
  run to tag mitochondrial contigs in ``{base}.final.fasta``
* `Flye <https://github.com/mikolmogorov/Flye>`_ -- long-read assembly (``pipeline_long``,
  ``pipeline_hybrid``)
* `Racon <https://github.com/lbcb-sci/racon>`_, `Polypolish <https://github.com/rrwick/Polypolish>`_,
  `pypolca <https://github.com/gbouras13/pypolca>`_ + freebayes -- polishing (``pipeline_long``,
  ``pipeline_hybrid``)
* Python packages: biopython, psutil (and matplotlib for ``depth`` plots)

**Only needed for optional steps or non-default options:**

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Tool
     - Used by
   * - `sourmash <https://sourmash.readthedocs.io/>`_ >= 4.8 + `branchwater plugin
       <https://github.com/sourmash-bio/sourmash_plugin_branchwater>`_
     - ``taxpurge`` (run by the pipelines only with ``-p/--phylum``)
   * - `fastp <https://github.com/OpenGene/fastp>`_, `Trimmomatic <https://github.com/usadellab/Trimmomatic>`_
     - ``trim --method fastp`` / ``trimmomatic``
   * - `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
     - ``filter --aligner bowtie2``
   * - `MEGAHIT <https://github.com/voutcn/megahit>`_, `Unicycler <https://github.com/rrwick/Unicycler>`_
     - ``assemble --method megahit`` / ``unicycler`` (also the pipelines' ``--method``)
   * - `NOVOPlasty <https://github.com/ndierckx/NOVOPlasty>`_
     - ``mito``
   * - `NCBI FCS-adaptor <https://github.com/ncbi/fcs>`_ (with singularity/apptainer or docker)
     - ``fcs_screen``
   * - `NCBI FCS-GX <https://github.com/ncbi/fcs-gx>`_
     - ``fcs_gx_purge`` (needs RAM for the whole ~470 GB database)
   * - `NextPolish2 <https://github.com/Nextomics/NextPolish2>`_ + yak
     - ``polish --method nextpolish2``
   * - pigz
     - faster read counting (falls back to gzip)

Checking the installation
=========================

``AAFTF dependency`` lists every external tool and Python package AAFTF uses, in two groups: those
the ``pipeline_*`` subcommands need with their default settings (BBTools and Java, SPAdes, BLAST+,
bwa, samtools, minimap2, Flye, Racon, Polypolish, pypolca, freebayes, mosdepth; biopython, psutil), and those only needed by optional steps
(``taxpurge``, ``polish``, ``mito``, ``fcs_screen``, ``fcs_gx_purge``) or non-default options
(e.g. ``filter --aligner bowtie2``, ``assemble --method megahit``). It exits with an error only
when something from the first group is missing; ``-q`` prints just the OK/MISSING status.

.. code-block:: bash

    AAFTF dependency

Setting up the reference database (``AAFTF_DB``)
=====================================================

Several subcommands (``filter``, ``vecscreen``, ``taxpurge``, ``fcs_screen``) need UniVec_Core, PhiX,
the ``contamdb`` contamination database built from them, sourmash taxonomy databases, and/or
the NCBI FCS-adaptor resources. Download them once with ``AAFTF database``; those subcommands then
read them from the database folder (and tell you what to download if something is missing). Without
``$AAFTF_DB`` this is ``~/.cache/aaftf``; because some databases are several GB, it is best to
point ``$AAFTF_DB`` at a folder with plenty of space (it can also list several folders separated
by ``:``, like ``$PATH``):

.. code-block:: bash

    export AAFTF_DB=/path/with/space/aaftf_db
    AAFTF database                                        # list databases
    AAFTF database required   # what the default pipeline needs (phix univec contamdb)

See :doc:`commands/database` for the list of databases and which subcommands use each (the
sourmash indices are large, so fetch only the one you need).

Runs launched inside the pixi-managed Docker/Singularity images default ``AAFTF_DB`` to
``/opt/aaftf_db``; bind-mount or ``-e AAFTF_DB=...`` to point at your own copy instead of
re-downloading it on every run.
