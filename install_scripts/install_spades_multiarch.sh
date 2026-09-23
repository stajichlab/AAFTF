#!/usr/bin/env bash
# Build SPAdes from source for two x86-64 micro-architecture levels and
# install per-command dispatch wrappers that pick the right build at runtime.
#
# Why: the bioconda SPAdes recipe patches spades_compile.sh to force
# CXXFLAGS=-march=x86-64-v3 (AVX2/BMI2/FMA). Every bioconda build since then
# dies with SIGILL (exit 132) on older CPUs that are only x86-64-v2, e.g.
# AMD Opteron "Abu Dhabi" (UCR HPCC c01-c30) and Intel Ivy Bridge (h01-h06).
# The crash happens at program start-up (a BMI2 `rorx` inside the bundled
# mimalloc), so no SPAdes command works on those nodes.
#
# What this script does:
#   1. Downloads the SPAdes source tarball and checks its sha256.
#   2. Builds it twice, into ${SPADES_ROOT}/v2 (-march=x86-64-v2) and
#      ${SPADES_ROOT}/v3 (-march=x86-64-v3). Each tree is self-contained:
#      spades_init.py locates bin/ and share/spades/ from its own realpath.
#   3. Checks that the v2 tree has no BMI2 instructions (objdump) and runs
#      `spades.py --test` from each tree the build host can execute.
#   4. Installs ${SPADES_ROOT}/spades-arch-exec, which tests /proc/cpuinfo
#      for the x86-64-v3 feature set and execs the matching tree.
#   5. Replaces every SPAdes command in ${BIN_DIR} (default $CONDA_PREFIX/bin,
#      which holds the bioconda v3-only binaries) with a two-line wrapper
#      around spades-arch-exec. AAFTF keeps calling `spades.py` unchanged.
#
# Override the runtime choice with SPADES_ARCH=v2 or SPADES_ARCH=v3.
#
# Build requirements: g++ >= 9.1, cmake, make, zlib and bzip2 headers
# (Debian/Ubuntu: build-essential cmake zlib1g-dev libbz2-dev), objdump
# (binutils), python3, curl. Run inside the activated pixi env.
#
# The compile runs with PATH=${SPADES_BUILD_PATH} (default: system dirs only).
# With the pixi env bin/ on PATH, CMake would also search the env prefix and
# could link SPAdes against conda libraries that are not on the runtime
# library path. `pixi run install-spades-multiarch` sets SPADES_BUILD_PATH to
# the full PATH because a pixi-only host may have cmake only in the env.

set -euo pipefail

SPADES_VERSION="${SPADES_VERSION:-4.3.0}"
SPADES_SHA256="${SPADES_SHA256:-09671ca39f9c6d2479d9fc168100bfd089b4a24002d51b815386d2b24d424456}"
SPADES_ROOT="${SPADES_ROOT:-/opt/spades}"
BIN_DIR="${BIN_DIR:-${CONDA_PREFIX:?CONDA_PREFIX is not set -- run inside a pixi/conda env or set BIN_DIR}/bin}"
NPROC="${NPROC:-${SLURM_CPUS_ON_NODE:-$(nproc)}}"
LEVELS=(v2 v3)
SPADES_BUILD_PATH="${SPADES_BUILD_PATH:-/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin}"

# SPAdes commands that get a dispatch wrapper (bioconda installs these into bin/).
COMMANDS=(
    spades.py metaspades.py plasmidspades.py metaplasmidspades.py rnaspades.py
    rnaviralspades.py metaviralspades.py coronaspades.py
    spades-hammer spades-ionhammer spades-core spades-corrector-core
    spades-gbuilder spades-gmapper spades-gsimplifier spades-kmercount
    spades-kmer-estimating spades-read-filter spades-convert-bin-to-fasta
    spades-gfa-split spades-bwa spades-hpc binspreader pathracer
    pathracer-seq-fs spaligner splitter
)

log() { printf '[install_spades_multiarch] %s\n' "$*"; }
die() { printf '[install_spades_multiarch] ERROR: %s\n' "$*" >&2; exit 1; }

for tool in cmake make g++ objdump strip python3 curl sha256sum; do
    PATH="${SPADES_BUILD_PATH}" command -v "$tool" >/dev/null 2>&1 \
        || die "required tool not found on SPADES_BUILD_PATH: $tool"
done

# Returns 0 when this host supports the x86-64-v3 feature set.
host_is_v3() {
    local flags f
    flags=" $(awk -F: '/^flags/ {print $2; exit}' /proc/cpuinfo) "
    for f in avx avx2 bmi1 bmi2 f16c fma abm movbe xsave; do
        [[ "$flags" == *" $f "* ]] || return 1
    done
}

BUILD_DIR=$(mktemp -d)
trap 'rm -rf "${BUILD_DIR}"' EXIT
cd "${BUILD_DIR}"

TARBALL="SPAdes-${SPADES_VERSION}.tar.gz"
log "downloading ${TARBALL}"
curl -fsSL -o "${TARBALL}" \
    "https://github.com/ablab/spades/releases/download/v${SPADES_VERSION}/${TARBALL}"
echo "${SPADES_SHA256}  ${TARBALL}" | sha256sum -c - >/dev/null \
    || die "sha256 mismatch for ${TARBALL}"
tar xzf "${TARBALL}"
SRC="${BUILD_DIR}/SPAdes-${SPADES_VERSION}"
[[ -d "${SRC}/src" ]] || die "unexpected tarball layout: ${SRC}/src missing"

for lvl in "${LEVELS[@]}"; do
    prefix="${SPADES_ROOT}/${lvl}"
    flags="-O3 -march=x86-64-${lvl} -mtune=generic -fcommon"
    log "building ${lvl} (${flags}) -> ${prefix} with ${NPROC} jobs"
    rm -rf "${prefix}"
    if ! PATH="${SPADES_BUILD_PATH}" CFLAGS="${flags}" \
        CXXFLAGS="${flags} -Wno-deprecated-declarations" \
        cmake -S "${SRC}/src" -B "${BUILD_DIR}/build_${lvl}" \
            -G "Unix Makefiles" \
            -DCMAKE_INSTALL_PREFIX="${prefix}" \
            -DSPADES_BUILD_INTERNAL=OFF \
            > "${BUILD_DIR}/cmake_${lvl}.log" 2>&1; then
        tail -50 "${BUILD_DIR}/cmake_${lvl}.log" >&2
        die "cmake configure failed for ${lvl}"
    fi
    if ! PATH="${SPADES_BUILD_PATH}" make -C "${BUILD_DIR}/build_${lvl}" -j "${NPROC}" install \
            > "${BUILD_DIR}/make_${lvl}.log" 2>&1; then
        tail -50 "${BUILD_DIR}/make_${lvl}.log" >&2
        die "make install failed for ${lvl}"
    fi
    strip "${prefix}"/bin/spades-* 2>/dev/null || true
    rm -rf "${BUILD_DIR}/build_${lvl}"
done

# The v2 tree must not contain BMI2 scalar instructions: those only appear
# when the whole translation unit is built for v3, and they are what raised
# SIGILL in the bioconda build. (Bundled zlib-ng keeps runtime-dispatched
# AVX2 kernels, so AVX2 vector opcodes alone are not a failure signal.)
log "checking v2 binaries for BMI2 instructions"
bad=0
for exe in "${SPADES_ROOT}"/v2/bin/*; do
    [[ -f "$exe" && -x "$exe" ]] || continue
    head -c 4 "$exe" | grep -q $'\x7fELF' || continue
    n=$(objdump -d --no-show-raw-insn "$exe" \
        | awk '$2 ~ /^(shlx|shrx|sarx|rorx|pdep|pext|mulx|bzhi)$/' | wc -l)
    if [[ "$n" -gt 0 ]]; then
        log "  $(basename "$exe"): ${n} BMI2 instructions"
        bad=1
    fi
done
[[ $bad -eq 0 ]] || die "v2 build contains x86-64-v3 code; check CFLAGS/CXXFLAGS"

for lvl in "${LEVELS[@]}"; do
    if [[ "$lvl" == v3 ]] && ! host_is_v3; then
        log "skipping v3 self-test: build host is not x86-64-v3"
        continue
    fi
    log "running spades.py --test for ${lvl}"
    # --test writes ./spades_test and rejects -o, so run it in its own dir.
    mkdir -p "${BUILD_DIR}/test_${lvl}"
    ( cd "${BUILD_DIR}/test_${lvl}" && "${SPADES_ROOT}/${lvl}/bin/spades.py" --test -t 2 \
        > "${BUILD_DIR}/test_${lvl}.log" 2>&1 ) \
        || { tail -30 "${BUILD_DIR}/test_${lvl}.log" >&2; die "spades.py --test failed for ${lvl}"; }
done

log "installing ${SPADES_ROOT}/spades-arch-exec"
cat > "${SPADES_ROOT}/spades-arch-exec" <<EOF
#!/usr/bin/env bash
# Run a SPAdes command from the build matching this CPU.
# Generated by AAFTF install_scripts/install_spades_multiarch.sh.
# Usage: spades-arch-exec <command> [args...]
# Set SPADES_ARCH=v2|v3 to force a build.
set -e
root="${SPADES_ROOT}"
cmd="\$1"; shift
lvl="\${SPADES_ARCH:-}"
if [[ -z "\$lvl" ]]; then
    lvl=v3
    flags=" \$(awk -F: '/^flags/ {print \$2; exit}' /proc/cpuinfo) "
    for f in avx avx2 bmi1 bmi2 f16c fma abm movbe xsave; do
        [[ "\$flags" == *" \$f "* ]] || { lvl=v2; break; }
    done
fi
exe="\$root/\$lvl/bin/\$cmd"
[[ -x "\$exe" ]] || { echo "spades-arch-exec: \$exe not found" >&2; exit 127; }
exec "\$exe" "\$@"
EOF
chmod 755 "${SPADES_ROOT}/spades-arch-exec"

log "installing wrappers into ${BIN_DIR}"
mkdir -p "${BIN_DIR}"
for cmd in "${COMMANDS[@]}"; do
    if [[ ! -x "${SPADES_ROOT}/v2/bin/${cmd}" || ! -x "${SPADES_ROOT}/v3/bin/${cmd}" ]]; then
        if [[ -e "${BIN_DIR}/${cmd}" ]]; then
            log "  WARNING: ${cmd} not built from source; leaving the conda (v3-only) copy"
        fi
        continue
    fi
    rm -f "${BIN_DIR}/${cmd}"
    printf '#!/bin/sh\nexec %s/spades-arch-exec %s "$@"\n' "${SPADES_ROOT}" "${cmd}" \
        > "${BIN_DIR}/${cmd}"
    chmod 755 "${BIN_DIR}/${cmd}"
done

touch "${SPADES_ROOT}/.multiarch_${SPADES_VERSION}"
log "done: $("${BIN_DIR}/spades.py" --version 2>&1) (this host uses $(host_is_v3 && echo v3 || echo v2))"
