"""Shared pytest fixtures for the AAFTF test suite."""

import gzip
import io
import textwrap
from argparse import Namespace
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# FASTA sequences with known properties
# seq1: 20 bp,  50 % GC
# seq2: 80 bp, 100 % GC
# seq3: 150 bp,   0 % GC
# Total: 250 bp, GC = (10+80+0)/250 = 36.0 %
# N50 = 150 (cumulative from largest reaches 50 % at first contig)
# N90 = 80  (cumulative reaches 90 % with first two contigs)
# ---------------------------------------------------------------------------
SEQ1 = "ATCGATCGATCGATCGATCG"   # 20 bp
SEQ2 = "GC" * 40                # 80 bp
SEQ3 = "AT" * 75                # 150 bp

SMALL_FASTA = f">seq1\n{SEQ1}\n>seq2\n{SEQ2}\n>seq3\n{SEQ3}\n"

# Telomere sequences (200 bp each).
# Default monomer: TAA[C]+  revcomp: [G]+TTA
_FWD  = "TAACCCTAACCC"           # matches TAA[C]+TAA[C]+
_REV  = "GGGTTAGGGTTA"           # matches [G]+TTA[G]+TTA
_BODY = "A" * 176                # neutral padding (176 bp)

# Forward telomere only (pattern in first 100 bp)
TELO_FWD_ONLY = _FWD + "A" * 188         # 200 bp
# Reverse telomere only (pattern in last 100 bp)
TELO_REV_ONLY = "A" * 188 + _REV        # 200 bp
# Telomere-to-telomere: both ends
TELO_T2T      = _FWD + _BODY + _REV     # 200 bp
# No telomere
TELO_NONE     = "AAATTTGGGCCC" * 16 + "AAAA" * 2   # 200 bp

TELOMERE_FASTA = (
    f">scaffold_t2t\n{TELO_T2T}\n"
    f">scaffold_fwd\n{TELO_FWD_ONLY}\n"
    f">scaffold_rev\n{TELO_REV_ONLY}\n"
    f">scaffold_none\n{TELO_NONE}\n"
)

# ---------------------------------------------------------------------------
# mosdepth summary data
# 9 nuclear contigs at 10x + 1 organellar contig at 1000x
# total_row mean ≈ 20.88x (length-weighted)
# SD of per-contig depths vs mean=20.88 ≈ 309.8
# threshold ≈ 950.3  →  scaffold_outlier (1000x) IS flagged
# ---------------------------------------------------------------------------
_NUCLEAR_ROWS = "\n".join(
    f"scaffold_{i}\t10000\t100000\t10.0\t0\t50" for i in range(1, 10)
)
MOSDEPTH_SUMMARY = (
    "chrom\tlength\tbases\tmean\tmin\tmax\n"
    + _NUCLEAR_ROWS
    + "\nscaffold_outlier\t1000\t1000000\t1000.0\t0\t2000"
    + "\ntotal\t91000\t1900000\t20.88\t0\t2000\n"
)

MOSDEPTH_DIST = (
    "total\t0\t1.00000\n"
    "total\t1\t0.98000\n"
    "total\t2\t0.96000\n"
    "total\t5\t0.90000\n"
)

# ---------------------------------------------------------------------------
# FASTQ helpers
# ---------------------------------------------------------------------------

def make_fastq_text(n_reads=10, read_len=100):
    """Return a FASTQ string with *n_reads* dummy records of *read_len* bp."""
    lines = []
    for i in range(n_reads):
        lines += [f"@read_{i}", "A" * read_len, "+", "I" * read_len]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# File fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fasta_file(tmp_path):
    """Small multi-sequence FASTA file with known properties."""
    p = tmp_path / "test.fasta"
    p.write_text(SMALL_FASTA)
    return p


@pytest.fixture
def gz_fasta_file(tmp_path):
    """Gzip-compressed version of the small FASTA file."""
    p = tmp_path / "test.fasta.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(SMALL_FASTA)
    return p


@pytest.fixture
def telomere_fasta_file(tmp_path):
    """FASTA file containing sequences with and without telomere repeats."""
    p = tmp_path / "telomere.fasta"
    p.write_text(TELOMERE_FASTA)
    return p


@pytest.fixture
def fastq_file(tmp_path):
    """Plain FASTQ file with 10 dummy reads."""
    p = tmp_path / "reads.fastq"
    p.write_text(make_fastq_text(10))
    return p


@pytest.fixture
def gz_fastq_file(tmp_path):
    """Gzip-compressed FASTQ file with 10 dummy reads."""
    p = tmp_path / "reads.fastq.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(make_fastq_text(10))
    return p


@pytest.fixture
def empty_file(tmp_path):
    """An empty file (zero bytes)."""
    p = tmp_path / "empty.txt"
    p.write_text("")
    return p


@pytest.fixture
def line_file(tmp_path):
    """A file with exactly 5 lines."""
    p = tmp_path / "lines.txt"
    p.write_text("line1\nline2\nline3\nline4\nline5\n")
    return p


@pytest.fixture
def mosdepth_summary_file(tmp_path):
    """mosdepth-format summary.txt with 9 nuclear contigs + 1 outlier."""
    p = tmp_path / "coverage.mosdepth.summary.txt"
    p.write_text(MOSDEPTH_SUMMARY)
    return p


@pytest.fixture
def mosdepth_dist_file(tmp_path):
    """mosdepth global distribution file."""
    p = tmp_path / "coverage.mosdepth.global.dist.txt"
    p.write_text(MOSDEPTH_DIST)
    return p


# ---------------------------------------------------------------------------
# args Namespace fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sort_args(tmp_path, fasta_file):
    """Minimal Namespace for sort.run()."""
    return Namespace(
        input=str(fasta_file),
        out=str(tmp_path / "sorted.fasta"),
        minlen=0,
        name="scaffold",
        debug=False,
        pipe=True,
    )


@pytest.fixture
def assess_args(tmp_path, fasta_file):
    """Minimal Namespace for assess.run()."""
    return Namespace(
        input=str(fasta_file),
        report=str(tmp_path / "report.txt"),
        telomere_monomer="TAA[C]+",
        telomere_n_repeat=2,
        debug=False,
        pipe=True,
    )
