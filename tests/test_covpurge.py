"""Unit tests for aaftf/covpurge.py (bwa/samtools are mocked)."""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.covpurge import run

pytestmark = pytest.mark.unit


def _ids(path):
    return [line[1:].strip() for line in Path(path).read_text().splitlines() if line.startswith(">")]


def _bedcov(covs):
    return [f"{k}\t0\t{length}\t{int(length * c)}\n" for k, (length, c) in covs.items()]


def _run(tmp_path, seqs, bedcov, **kw):
    """Run covpurge with ``samtools bedcov`` returning ``bedcov``; return (kept ids, run_cmd, execute, align mocks)."""
    asm = tmp_path / "asm.fasta"
    asm.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    out = tmp_path / "out.fasta"
    with (
        patch("aaftf.covpurge.run_cmd") as rc,
        patch("aaftf.covpurge.execute", side_effect=lambda *a, **k: iter(bedcov)) as ex,
        patch("aaftf.covpurge.align_to_sorted_bam") as bam,
    ):
        kw.setdefault("read1", str(tmp_path / "r1.fq"))
        run(input=str(asm), outfile=str(out), workdir=str(tmp_path / "wd"), **kw)
    return _ids(out), rc, ex, bam


@pytest.mark.parametrize("paired", [False, True])
def test_coverage_filter(tmp_path, paired):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    seqs = {"big": "A" * 1000, "mid": "C" * 500, "low": "G" * 400, "edge": "T" * 300}
    # N50 = 500 bp, so N50 contigs are big+mid (avg 100X): cutoff 5% -> 5X; <= cutoff dropped
    bed = _bedcov({"big": (1000, 100), "mid": (500, 100), "low": (400, 5.5), "edge": (300, 5)})
    kw = {"read1": str(r1), "read2": str(r2) if paired else None}
    kept, rc, ex, bam = _run(tmp_path, seqs, bed, **kw)
    assert kept == ["big", "mid", "low"]
    assert rc.call_args.args[0] == ["bwa", "index", "assembly.fasta"]
    bwa = bam.call_args.args[0]
    assert bwa[:5] == ["bwa", "mem", "-t", "1", "assembly.fasta"]
    assert bwa[5:] == ([str(r1.resolve()), str(r2.resolve())] if paired else [str(r1.resolve())])
    assert ex.call_args.args[0] == ["samtools", "bedcov", "assembly.bed", "remapped.bam"]


def test_mincovpct_changes_cutoff(tmp_path):
    seqs = {"big": "A" * 1000, "mid": "C" * 500}
    kept, *_ = _run(tmp_path, seqs, _bedcov({"big": (1000, 100), "mid": (500, 20)}), mincovpct=20)
    assert kept == ["big"]


def test_no_n50_coverage_warns_and_skips(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    seqs = {"big": "A" * 1000, "low": "G" * 400}
    kept, *_ = _run(tmp_path, seqs, _bedcov({"low": (400, 0)}))
    assert kept == ["big", "low"]
    assert "skipping the low-coverage filter" in caplog.text


def test_existing_bam_skips_alignment(tmp_path):
    (tmp_path / "wd").mkdir()
    (tmp_path / "wd" / "remapped.bam").write_text("")
    kept, rc, _, bam = _run(tmp_path, {"contig": "A" * 100}, _bedcov({"contig": (100, 10)}))
    rc.assert_not_called()
    bam.assert_not_called()
    assert kept == ["contig"]


def test_pipe_hides_hint(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    seqs, bed = {"contig": "A" * 100}, _bedcov({"contig": (100, 10)})
    _run(tmp_path, seqs, bed, pipe=True)
    assert "AAFTF rmdup" not in caplog.text
    _run(tmp_path, seqs, bed, pipe=False)
    assert "AAFTF rmdup" in caplog.text
