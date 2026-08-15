#!/usr/bin/env bash
# Build bowtie2 from source to fix AVX2 runtime fallback issues in the conda
# binary. Mirrors the approach used in funannotate.
#
# The conda-packaged bowtie2 2.5.5 fails to launch the x86-64-v3 (AVX2)
# optimized code at runtime, falling back to slower baseline x86-64
# performance (warns: "Failed to launch x86-64-v3 version, staying with
# default"). Building from source with the local toolchain resolves this.
#
# This script:
# 1. Clones bowtie2 v2.5.5 from GitHub
# 2. Compiles with the system toolchain (g++, make)
# 3. Installs binaries to $CONDA_PREFIX/bin, shadowing the conda package
# 4. Is idempotent: safe to run multiple times
#
# Idempotency uses a marker file placed after a successful source build. A
# runtime check alone is NOT reliable because the conda binary's SIMD runtime
# dispatch emits the x86-64-v3 fallback warning non-deterministically, so it
# could falsely conclude the fix was already applied.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set -- run inside a pixi/conda env}"

BT2_VERSION="2.5.5"
BT2_BIN_DIR="${CONDA_PREFIX}/bin"
BT2_MARKER="${BT2_BIN_DIR}/.bowtie2_sourcebuild"

# If we previously built from source, confirm the install is still clean, then
# exit. Otherwise (first run on a fresh conda env) always rebuild.
if [ -f "${BT2_MARKER}" ]; then
    echo "[pixi_install_bowtie2] Source build already installed, checking for runtime issues..."
    if ! "${BT2_BIN_DIR}/bowtie2" --version 2>&1 | grep -q "WARNING"; then
        echo "[pixi_install_bowtie2] No x86-64-v3 warnings detected, installation OK"
        exit 0
    fi
    echo "[pixi_install_bowtie2] Marker present but AVX2 warnings detected, rebuilding..."
else
    echo "[pixi_install_bowtie2] Building bowtie2 ${BT2_VERSION} from source..."
fi

# Create temporary build directory
BT2_BUILD_DIR=$(mktemp -d)
trap "rm -rf '${BT2_BUILD_DIR}'" EXIT

cd "${BT2_BUILD_DIR}"

echo "[pixi_install_bowtie2] Cloning bowtie2 v${BT2_VERSION}..."
git clone --depth 1 --branch v${BT2_VERSION} \
    https://github.com/BenLangmead/bowtie2.git bowtie2

cd bowtie2

echo "[pixi_install_bowtie2] Compiling with local toolchain..."
# Use make with parallel jobs (use SLURM_CPUS_ON_NODE if available, otherwise 16).
# Requires a C++ compiler that understands -march=x86-64-v3 (GCC >= 8) and the
# "x86-64-v3" __builtin_cpu_supports argument (GCC >= 12) -- e.g. the GCC 12+
# shipped by Debian bookworm (/usr/bin/g++) via build-essential.
NPROC=${SLURM_CPUS_ON_NODE:-16}
if ! make -j "${NPROC}" > /tmp/bowtie2_build.log 2>&1; then
    echo "[pixi_install_bowtie2] ERROR: Build failed. Last 50 lines of build log:" >&2
    tail -50 /tmp/bowtie2_build.log >&2
    exit 1
fi

echo "[pixi_install_bowtie2] Build successful, installing binaries..."

# Install everything into $CONDA_PREFIX/bin (PREFIX -> $CONDA_PREFIX, bindir ->
# $CONDA_PREFIX/bin). This copies the launcher scripts plus bowtie2-align-*,
# bowtie2-build-*, bowtie2-inspect-* AND the -v256 (x86-64-v3/AVX2) variants.
# The conda package omits the -v256 binaries, which is exactly why its runtime
# dispatch fails with "Failed to launch x86-64-v3 version, staying with
# default" and silently drops to baseline SSE2. Shipping them restores AVX2.
if ! make install PREFIX="${CONDA_PREFIX}" > /tmp/bowtie2_install.log 2>&1; then
    echo "[pixi_install_bowtie2] ERROR: make install failed. Last 50 lines of install log:" >&2
    tail -50 /tmp/bowtie2_install.log >&2
    exit 1
fi

echo "[pixi_install_bowtie2] Installation complete. Verifying..."

# Mark this env as source-built (used by the idempotency check above)
touch "${BT2_MARKER}"

# Test that new binary works and check for warnings
BT2_OUTPUT=$("${BT2_BIN_DIR}/bowtie2" --version 2>&1)
echo "[pixi_install_bowtie2] Version: $BT2_OUTPUT"

if echo "$BT2_OUTPUT" | grep -q "WARNING"; then
    echo "[pixi_install_bowtie2] WARNING: x86-64-v3 fallback still occurring (may be expected)" >&2
else
    echo "[pixi_install_bowtie2] SUCCESS: No AVX2 fallback warnings detected"
fi

exit 0
