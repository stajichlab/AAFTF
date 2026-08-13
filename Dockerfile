# =============================================================================
# AAFTF - Automatic Assembly For The Fungi
# Docker image built with pixi-managed conda + PyPI dependencies.
#
# Build (from source, dev/editable install). AAFTF_VERSION must be supplied
# explicitly: .dockerignore excludes .git from the build context, so
# setuptools-scm cannot derive the version inside the image the way it does
# in a real source checkout (see AAFTF/__version__.py) — pass the host's
# git-derived version through instead:
#   docker build --build-arg AAFTF_VERSION=$(git describe --tags --always) \
#       -t aaftf:latest .
#
# Build for a specific tagged release (uses the "release" pixi environment):
#   docker build --build-arg PIXI_ENV=release --build-arg AAFTF_VERSION=0.6.2 \
#       -t aaftf:v0.6.2 .
#
# Run:
#   docker run --rm aaftf:latest AAFTF --help
#
# Run with data and external DB:
#   docker run --rm \
#     -v /path/to/data:/data \
#     -v /path/to/aaftf_db:/opt/aaftf_db \
#     aaftf:latest AAFTF trim --left /data/R1.fq.gz --right /data/R2.fq.gz
# =============================================================================

FROM debian:bookworm-slim

# Which pixi environment to activate (default = editable local install).
# Pass --build-arg PIXI_ENV=release to install from the pinned git tag instead.
ARG PIXI_ENV=default

# setuptools-scm can't see git history inside the build context (.git is
# excluded via .dockerignore), so the version must be supplied explicitly.
ARG AAFTF_VERSION=0.0.0+unknown
ENV SETUPTOOLS_SCM_PRETEND_VERSION_FOR_AAFTF=${AAFTF_VERSION}

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
        locales-all && \
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
#    pixi.toml + pixi.lock are copied first so dependency installation is
#    cached independently of source-code changes.
# ---------------------------------------------------------------------------
WORKDIR /opt/AAFTF
COPY pixi.toml pixi.lock ./

# Copy rest of the source (needed before pixi install because the default
# environment installs AAFTF as an editable PyPI package from ".")
COPY . .

# ---------------------------------------------------------------------------
# 4. Install conda + PyPI dependencies via pixi (locked / reproducible)
# ---------------------------------------------------------------------------
RUN pixi install --environment "${PIXI_ENV}" --locked && \
    # Export activation env-vars so every subsequent CMD can use them.
    pixi shell-hook --environment "${PIXI_ENV}" --shell bash \
        | grep '^export ' > /opt/aaftf_activate.sh && \
    chmod +x /opt/aaftf_activate.sh && \
    pixi clean cache --yes

# ---------------------------------------------------------------------------
# 5. Apply the patched polca.sh (compatible with newer samtools)
# ---------------------------------------------------------------------------
RUN source /opt/aaftf_activate.sh && \
    POLCA_INSTALLED=$(find /opt/AAFTF/.pixi/envs -name "polca.sh" | head -1) && \
    if [ -n "$POLCA_INSTALLED" ]; then \
        echo "Replacing $POLCA_INSTALLED with patched version"; \
        cp /opt/AAFTF/patches/polca.sh "$POLCA_INSTALLED"; \
        chmod +x "$POLCA_INSTALLED"; \
    else \
        echo "WARNING: polca.sh not found — skipping polca patch"; \
    fi

# ---------------------------------------------------------------------------
# 6. Smoke test
# ---------------------------------------------------------------------------
RUN source /opt/aaftf_activate.sh && AAFTF --version

# ---------------------------------------------------------------------------
# 7. Cleanup build artefacts to reduce image size
#    Keep /opt/pixi intact — pixi manages the conda env and removing it can
#    break activation.  Only purge the download cache.
# ---------------------------------------------------------------------------
RUN rm -rf /root/.cache

# ---------------------------------------------------------------------------
# 8. Entrypoint: source the activation script then exec the user command
# ---------------------------------------------------------------------------
RUN printf '#!/bin/bash\nset -e\nsource /opt/aaftf_activate.sh\nexec "$@"\n' \
        > /usr/local/bin/docker-entrypoint.sh && \
    chmod +x /usr/local/bin/docker-entrypoint.sh

# Default AAFTF_DB location (override at runtime with -e or -v).
ENV AAFTF_DB="/opt/aaftf_db"

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["AAFTF", "--help"]
