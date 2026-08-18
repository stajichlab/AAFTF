#!/usr/bin/env bash
# =============================================================================
# build_sif.sh — Build the AAFTF Singularity/Apptainer SIF from AAFTF.def
#
# Tested on UCR HPCC with apptainer 1.4.5 doing an UNPRIVILEGED build
# (user namespaces; no root / sudo / --fakeroot needed). This is the only
# non-root path that works here: singularity-ce 4.3.2 refuses unprivileged
# def builds ("--remote, --fakeroot, or the proot command are required").
#
# The produced SIF is fully self-contained: pixi env (conda + PyPI deps),
# AAFTF editable install, the bowtie2 source build (AVX2 / -v256 fix), and
# /etc/profile.d/aaftf.sh so AAFTF is on PATH even in login shells.
#
# Usage:
#   ./build_sif.sh                       # full build, bakes in sourmash DB
#   ./build_sif.sh --skip-db             # no sourmash DB (CI-style, faster)
#   ./build_sif.sh --out /path/to.sif --skip-db
#   ./build_sif.sh --git-ref v0.7.0      # build a specific tag/release
#
# Run (with singularity OR apptainer):
#   singularity exec AAFTF.sif AAFTF trim --left R1.fq.gz --right R2.fq.gz
#   singularity run  --bind /data:/data AAFTF.sif  AAFTF --help
#   singularity exec AAFTF.sif bash -lc 'AAFTF --version'   # login-shell test
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="${SCRIPT_DIR}/AAFTF.def"

OUT="${SCRIPT_DIR}/AAFTF.sif"
GIT_REF="main"
SKIP_DB=0
BUILD_ARGS=()

log() { printf '[build_sif] %s\n' "$*"; }
die() { printf '[build_sif] ERROR: %s\n' "$*" >&2; exit 1; }

usage() {
    sed -n '2,30p' "${BASH_SOURCE[0]}"
    cat <<EOF

Options:
  -o, --out PATH       Output .sif path                    (default: $OUT)
  -r, --git-ref REF    Git branch/tag of AAFTF to build    (default: $GIT_REF)
      --skip-db        Do not bake the sourmash DB into the image
      --fakeroot       Try a --fakeroot build (needs subuid mapping)
  -h, --help           Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--out)        OUT="$2"; shift 2 ;;
        -r|--git-ref)    GIT_REF="$2"; shift 2 ;;
        --skip-db)       SKIP_DB=1; shift ;;
        --fakeroot)      BUILD_ARGS+=(--fakeroot); shift ;;
        -h|--help)       usage; exit 0 ;;
        *) die "unknown argument: $1 (try --help)" ;;
    esac
done

[[ -f "$DEF_FILE" ]] || die "definition file not found: $DEF_FILE"

# -----------------------------------------------------------------------------
# Locate apptainer (module, absolute install path, or PATH)
# -----------------------------------------------------------------------------
APPTAINER_BIN=""
for cand in \
    "$(command -v apptainer 2>/dev/null || true)" \
    /usr/local/bin/apptainer \
    /opt/linux/rocky/8.x/x86_64/pkgs/apptainer/1.4.5/bin/apptainer; do
    if [[ -x "$cand" ]]; then APPTAINER_BIN="$cand"; break; fi
done

if [[ -z "$APPTAINER_BIN" ]]; then
    if command -v module >/dev/null 2>&1 && module load apptainer/1.4.5 2>/dev/null \
            && command -v apptainer >/dev/null 2>&1; then
        APPTAINER_BIN="$(command -v apptainer)"
    fi
fi
[[ -z "$APPTAINER_BIN" ]] && die "apptainer not found (try: module load apptainer/1.4.5)"
log "using apptainer: $APPTAINER_BIN ($("$APPTAINER_BIN" --version))"

[[ $SKIP_DB -eq 1 ]] && BUILD_ARGS+=(--build-arg skip_db_download=1)
BUILD_ARGS+=(--build-arg "git_ref=${GIT_REF}")

if [[ $SKIP_DB -eq 1 ]]; then
    log "build from ref '${GIT_REF}' -> ${OUT} (skip sourmash DB)"
else
    log "build from ref '${GIT_REF}' -> ${OUT} (WITH sourmash DB download)"
fi

rm -f "$OUT"
set -x
"$APPTAINER_BIN" build "${BUILD_ARGS[@]}" "$OUT" "$DEF_FILE"
set +x

log "done: $OUT"
if [[ -f "$OUT" ]]; then
    log "verify:"
    log "  apptainer exec $OUT bowtie2 --version"
    log "  apptainer exec $OUT bash -lc 'AAFTF --version'"
else
    die "build did not produce $OUT"
fi
