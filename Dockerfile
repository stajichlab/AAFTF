# =============================================================================
# AAFTF - Automatic Assembly For The Fungi
# Docker image built with pixi-managed conda + PyPI dependencies.
#
# Build (from source, dev/editable install). AAFTF_VERSION must be supplied
# explicitly: .dockerignore excludes .git from the build context, so hatch-vcs
# cannot derive the version inside the image. The --build-arg (a tag like
# v0.7.0 or `git describe` output like v0.7.0-beta.4-63-g7c3c423, converted to
# PEP 440 in step 4) is given to hatch-vcs as SETUPTOOLS_SCM_PRETEND_VERSION, so
# it ends up in the installed package metadata that `AAFTF --version` reads:
#   docker build --build-arg AAFTF_VERSION=$(git describe --tags --always) \
#       -t aaftf:latest .
#
# Without this arg the image will report v0.0.0+unknown.
#
# The image uses the pixi "complete" environment (every tool AAFTF can use). For a smaller
# image with only what the default pipeline needs:
#   docker build --build-arg PIXI_ENV=default --build-arg AAFTF_VERSION=... -t aaftf:slim .
#
# Run:
#   docker run --rm aaftf:latest AAFTF --help
#
# Run with data and external DB:
#   docker run --rm \
#     -v /path/to/data:/data \
#     -v /path/to/aaftf_db:/opt/aaftf_db \
#     aaftf:latest AAFTF trim --read1 /data/R1.fq.gz --read2 /data/R2.fq.gz
# =============================================================================

FROM ubuntu:noble

LABEL org.opencontainers.image.description="Automatic Assembly For The Fungi"

# Which pixi environment to install: "complete" (default-pipeline tools + every optional
# tool) or "default" (only what `AAFTF pipeline_short` needs with its default settings).
ARG PIXI_ENV=complete

# hatch-vcs can't see git history inside the build context (.git is excluded
# via .dockerignore), so the version must be supplied explicitly; it is used
# during `pixi install` (step 4) when aaftf is installed.
ARG AAFTF_VERSION=0.0.0+unknown

LABEL org.opencontainers.image.title="AAFTF"
LABEL org.opencontainers.image.description="Automatic Assembly For The Fungi"
LABEL org.opencontainers.image.url="https://github.com/stajichlab/AAFTF"
LABEL org.opencontainers.image.source="https://github.com/stajichlab/AAFTF"
LABEL org.opencontainers.image.licenses="MIT"

# Use bash for all RUN steps so we can source shell scripts.
SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

# ---------------------------------------------------------------------------
# 1. System dependencies
# ---------------------------------------------------------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        wget \
        bzip2 \
        perl \
        libgomp1 \
        procps \
        git \
        tar \
        pigz \
        aria2 \
        rsync \
        locales \
        locales-all \
        build-essential \
        cmake \
        libbz2-dev \
        binutils \
        patch \
        python3 \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Set locale for bioinformatics tools that require it (avoids
# "bash: warning: setlocale: LC_ALL: cannot change locale" when the host
# exports LC_ALL=en_US.UTF-8 but no en_US locale is generated in the image).
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    LANGUAGE=en_US:en

# ---------------------------------------------------------------------------
# 2. Install pixi
# ---------------------------------------------------------------------------
ENV PIXI_HOME=/opt/pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash
ENV PATH="${PIXI_HOME}/bin:${PATH}"

# ---------------------------------------------------------------------------
# 3. Copy repository into the image
#    pyproject.toml (pixi config lives under [tool.pixi.*]) + pixi.lock are
#    copied first so dependency installation is cached independently of
#    source-code changes.
# ---------------------------------------------------------------------------
WORKDIR /opt/AAFTF
COPY pyproject.toml pixi.lock ./

# Copy rest of the source (needed before pixi install because the default
# environment installs AAFTF as an editable PyPI package from ".")
COPY . .

# ---------------------------------------------------------------------------
# 4. Install conda + PyPI dependencies via pixi (locked / reproducible)
#    hatch-vcs only honours the generic SETUPTOOLS_SCM_PRETEND_VERSION (not the
#    _FOR_AAFTF form) and needs a PEP 440 version, so drop a leading "v" and
#    turn `git describe`'s "-g<hash>" into a "+g<hash>" local version.
# ---------------------------------------------------------------------------
RUN ver="${AAFTF_VERSION#v}" && \
    export SETUPTOOLS_SCM_PRETEND_VERSION="${ver/-g/+g}" && \
    pixi install --environment "${PIXI_ENV}" --locked && \
    # Export activation env-vars so every subsequent CMD can use them.
    pixi shell-hook --environment "${PIXI_ENV}" --shell bash \
        | grep '^export ' > /opt/aaftf_activate.sh && \
    chmod +x /opt/aaftf_activate.sh && \
    source /opt/aaftf_activate.sh && \
    python3 -c "import sys; from packaging.version import Version; from aaftf import __version__ as v; sys.exit(0 if Version(v) == Version(sys.argv[1]) else f'AAFTF reports {v}, expected {sys.argv[1]}')" "${SETUPTOOLS_SCM_PRETEND_VERSION}" && \
    pixi clean cache --yes

# ---------------------------------------------------------------------------
# 5. Build bowtie2 from source (fixes the conda binary's AVX2/x86-64-v3
#    runtime fallback). The conda bowtie2 2.5.5 fails to launch its v3 build
#    ("Failed to launch x86-64-v3 version, staying with default"), silently
#    dropping to baseline SSE2. install_scripts/pixi_install_bowtie2.sh rebuilds
#    it with the container toolchain and `make install PREFIX=$CONDA_PREFIX`,
#    which also ships the -v256 (AVX2) variants the conda package omits.
#    Mirrors the funannotate Dockerfile approach, and keeps pixi/Docker builds
#    from drifting apart.
# ---------------------------------------------------------------------------
RUN source /opt/aaftf_activate.sh && \
    bash install_scripts/pixi_install_bowtie2.sh && \
    test -x "${CONDA_PREFIX}/bin/bowtie2" && \
    test -x "${CONDA_PREFIX}/bin/bowtie2-align-s-v256"

# ---------------------------------------------------------------------------
# 6. Build SPAdes for x86-64-v2 AND x86-64-v3 with runtime dispatch.
#    The bioconda SPAdes recipe forces -march=x86-64-v3, so its binaries die
#    with SIGILL (exit 132) on older CPUs (e.g. AMD Opteron / Intel Ivy
#    Bridge HPC nodes). install_spades_multiarch.sh builds both levels under
#    /opt/spades/{v2,v3} and replaces the conda SPAdes commands in the env
#    bin/ with wrappers that pick the build matching the CPU at run time.
#    Force a build with SPADES_ARCH=v2|v3.
# ---------------------------------------------------------------------------
RUN source /opt/aaftf_activate.sh && \
    bash install_scripts/install_spades_multiarch.sh && \
    test -x /opt/spades/v2/bin/spades-hammer && \
    test -x /opt/spades/v3/bin/spades-hammer

# ---------------------------------------------------------------------------
# 7. Rebuild MEGAHIT from source without the bioconda -march=x86-64-v3 flag.
#    MEGAHIT picks one of three core binaries at run time (BMI2 / POPCNT /
#    baseline); the bioconda build puts v3 code into all three, so they die
#    with SIGILL on older CPUs. See install_scripts/install_megahit_from_source.sh.
# ---------------------------------------------------------------------------
RUN source /opt/aaftf_activate.sh && \
    bash install_scripts/install_megahit_from_source.sh && \
    test -f "${CONDA_PREFIX}/bin/.megahit_sourcebuild_1.2.9"

# ---------------------------------------------------------------------------
# 8. Smoke test
# ---------------------------------------------------------------------------
RUN source /opt/aaftf_activate.sh && AAFTF --version && \
    command -v aaftf_polca.sh fix_consensus_from_vcf.pl

# ---------------------------------------------------------------------------
# 9. Cleanup build artefacts to reduce image size
#    Keep /opt/pixi intact — pixi manages the conda env and removing it can
#    break activation.  Only purge the download cache.
# ---------------------------------------------------------------------------
RUN rm -rf /root/.cache

# ---------------------------------------------------------------------------
# 10. Entrypoint: source the activation script then exec the user command
# ---------------------------------------------------------------------------
RUN printf '#!/bin/bash\nset -e\nsource /opt/aaftf_activate.sh\nexec "$@"\n' \
        > /usr/local/bin/docker-entrypoint.sh && \
    chmod +x /usr/local/bin/docker-entrypoint.sh

# ---------------------------------------------------------------------------
# 11. Bake the pixi env bin dir onto PATH (mirrors funannotate's
#     `ENV PATH="/venv/bin:..."`). Docker execution gets this via the
#     ENTRYPOINT sourcing aaftf_activate.sh, but Singularity/Apptainer SIFs
#     converted from this image do NOT run the Docker ENTRYPOINT, and login
#     shells (SLURM/Nextflow use `bash -l`) reset PATH via /etc/profile.
#     Hardcode the same value aaftf_activate.sh exports so bowtie2/AAFTF are
#     found in all execution modes.
# ---------------------------------------------------------------------------
ENV PATH="/opt/AAFTF/.pixi/envs/${PIXI_ENV}/bin:/opt/pixi/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Default AAFTF_DB location (override at runtime with -e or -v).
ENV AAFTF_DB="/opt/aaftf_db"

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["AAFTF", "--help"]
