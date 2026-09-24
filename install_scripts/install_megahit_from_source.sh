#!/usr/bin/env bash
# Rebuild MEGAHIT from source so its own CPU dispatch works on x86-64-v2 CPUs.
#
# Why: MEGAHIT ships three core binaries and its `megahit` wrapper picks one at
# run time (`megahit_core checkcpu` / `checkpopcnt`):
#   megahit_core             -mbmi2 -mpopcnt   (BMI2 CPUs)
#   megahit_core_popcnt      -mpopcnt
#   megahit_core_no_hw_accel (baseline)
# The bioconda recipe adds a global CXXFLAGS=-march=x86-64-v3, so all three
# cores contain AVX2/BMI2 code and die with SIGILL on older CPUs, e.g. AMD
# Opteron "Abu Dhabi" (UCR HPCC c01-c30). Building without that flag restores
# upstream's dispatch: the popcnt / no_hw_accel cores then run on those nodes.
#
# What this script does:
#   1. Downloads MEGAHIT ${MEGAHIT_VERSION} and checks its sha256.
#   2. Applies the patches bioconda uses (patches/megahit/*.patch).
#   3. Builds with -O3 -mtune=generic (no -march) and installs into
#      ${MEGAHIT_PREFIX} (default $CONDA_PREFIX), replacing the conda copies of
#      megahit, megahit_toolkit and the three megahit_core* binaries.
#   4. Fails if megahit_core_popcnt or megahit_core_no_hw_accel contain BMI2
#      or AVX2 instructions, and runs `megahit --test`.
#
# Build requirements: g++, cmake, make, zlib headers (Debian/Ubuntu:
# build-essential cmake zlib1g-dev), objdump (binutils), patch, python3, curl.
# The compile runs with PATH=${MEGAHIT_BUILD_PATH} (default: system dirs only),
# so CMake does not pick up conda libraries from the pixi env prefix.

set -euo pipefail

MEGAHIT_VERSION="${MEGAHIT_VERSION:-1.2.9}"
MEGAHIT_SHA256="${MEGAHIT_SHA256:-09026eb07cc4e2d24f58b0a13f7a826ae8bb73da735a47cb1cbe6e4693118852}"
MEGAHIT_PREFIX="${MEGAHIT_PREFIX:-${CONDA_PREFIX:?CONDA_PREFIX is not set -- run inside a pixi/conda env or set MEGAHIT_PREFIX}}"
MEGAHIT_BUILD_PATH="${MEGAHIT_BUILD_PATH:-/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin}"
PATCH_DIR="${PATCH_DIR:-$(pwd)/patches/megahit}"
NPROC="${NPROC:-${SLURM_CPUS_ON_NODE:-$(nproc)}}"

log() { printf '[install_megahit_from_source] %s\n' "$*"; }
die() { printf '[install_megahit_from_source] ERROR: %s\n' "$*" >&2; exit 1; }

for tool in cmake make g++ objdump patch python3 curl sha256sum; do
    PATH="${MEGAHIT_BUILD_PATH}" command -v "$tool" >/dev/null 2>&1 \
        || die "required tool not found on MEGAHIT_BUILD_PATH: $tool"
done
compgen -G "${PATCH_DIR}/*.patch" >/dev/null || die "no patches found in ${PATCH_DIR} (run from the AAFTF repo root or set PATCH_DIR)"

BUILD_DIR=$(mktemp -d)
trap 'rm -rf "${BUILD_DIR}"' EXIT
cd "${BUILD_DIR}"

TARBALL="megahit-${MEGAHIT_VERSION}.tar.gz"
log "downloading ${TARBALL}"
curl -fsSL -o "${TARBALL}" "https://github.com/voutcn/megahit/archive/v${MEGAHIT_VERSION}.tar.gz"
echo "${MEGAHIT_SHA256}  ${TARBALL}" | sha256sum -c - >/dev/null \
    || die "sha256 mismatch for ${TARBALL}"
tar xzf "${TARBALL}"
SRC="${BUILD_DIR}/megahit-${MEGAHIT_VERSION}"

for p in "${PATCH_DIR}"/*.patch; do
    log "applying $(basename "$p")"
    patch -d "${SRC}" -p1 -s < "$p" || die "patch failed: $p"
done

flags="-O3 -mtune=generic -std=c++14 -Wno-deprecated-declarations -Wno-unused-but-set-variable"
log "building (${flags}) -> ${MEGAHIT_PREFIX} with ${NPROC} jobs"
if ! PATH="${MEGAHIT_BUILD_PATH}" CXXFLAGS="${flags}" \
    cmake -S "${SRC}" -B "${BUILD_DIR}/build" -G "Unix Makefiles" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${MEGAHIT_PREFIX}" \
        > "${BUILD_DIR}/cmake.log" 2>&1; then
    tail -50 "${BUILD_DIR}/cmake.log" >&2
    die "cmake configure failed"
fi
if ! PATH="${MEGAHIT_BUILD_PATH}" make -C "${BUILD_DIR}/build" -j "${NPROC}" install \
        > "${BUILD_DIR}/make.log" 2>&1; then
    tail -50 "${BUILD_DIR}/make.log" >&2
    die "make install failed"
fi
PATH="${MEGAHIT_BUILD_PATH}" strip "${MEGAHIT_PREFIX}"/bin/megahit_core* 2>/dev/null || true

# The fallback cores must be free of BMI2 and AVX2 code, otherwise the
# wrapper's dispatch cannot help on older CPUs.
log "checking fallback cores for BMI2/AVX2 instructions"
for core in megahit_core_popcnt megahit_core_no_hw_accel; do
    exe="${MEGAHIT_PREFIX}/bin/${core}"
    [[ -x "$exe" ]] || die "${exe} was not installed"
    n=$(objdump -d --no-show-raw-insn "$exe" \
        | awk '$2 ~ /^(shlx|shrx|sarx|rorx|pdep|pext|mulx|bzhi|vpbroadcast[bwdq]|vinserti128|vperm2i128|vpermq|vpgather[dq][dq])$/' | wc -l)
    [[ "$n" -eq 0 ]] || die "${core} contains ${n} BMI2/AVX2 instructions"
done

log "running megahit --test"
( cd "${BUILD_DIR}" && "${MEGAHIT_PREFIX}/bin/megahit" --test -t 2 -o "${BUILD_DIR}/megahit_test" \
    > "${BUILD_DIR}/test.log" 2>&1 ) \
    || { tail -30 "${BUILD_DIR}/test.log" >&2; die "megahit --test failed"; }

touch "${MEGAHIT_PREFIX}/bin/.megahit_sourcebuild_${MEGAHIT_VERSION}"
log "done: $("${MEGAHIT_PREFIX}/bin/megahit" --version 2>&1)"
