"""Unit tests for aaftf/vecscreen.py (blastn is mocked)."""

import logging
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf import vecscreen

pytestmark = pytest.mark.unit


def _read_fasta(path):
    """Return {id: seq} for a FASTA file."""
    recs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = ""
        elif name:
            recs[name] += line.strip()
    return recs


def _vec_row(qid, qstart, qend, score, qlen=1000):
    """A row of the vecscreen report (``_OUTFMT`` columns)."""
    return "\t".join([qid, "vec", "100", str(abs(qend - qstart) + 1), str(qstart), str(qend), "1e-5", "40", str(score), str(qlen)])


# --- _read_vector_hits --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("qstart", "qend", "score", "kept"),
    [
        (1, 30, 16, True),  # at the 5' tip: >= 16
        (1, 30, 15, False),
        (971, 1000, 16, True),  # at the 3' tip
        (10, 40, 19, True),  # terminal (within 50 bp), not at the tip: >= 19
        (10, 40, 18, False),
        (950, 960, 19, True),  # 3' terminal: 1000 - 960 < 50
        (500, 530, 30, True),  # internal: Strong only
        (500, 530, 29, False),  # internal moderate no longer splits
        (51, 80, 29, False),  # 51 > 50 -> internal
        (600, 500, 30, True),  # minus-strand coordinates are sorted
    ],
)
def test_read_vector_hits_thresholds(tmp_path, qstart, qend, score, kept):
    tab = tmp_path / "r.tab"
    tab.write_text(_vec_row("c1", qstart, qend, score) + "\n")
    hits = vecscreen._read_vector_hits(str(tab))
    assert hits == ({"c1": [tuple(sorted((qstart, qend)))]} if kept else {})


# --- _vector_regions ----------------------------------------------------------------------------


def test_vector_regions_merge_and_extend_to_ends():
    # 5' hit within 50 bp is extended to the start; hits < 50 bp apart merge; 3' hit extended to the end
    assert vecscreen._vector_regions([(20, 40), (500, 520), (560, 580), (300, 320), (310, 330), (960, 975)], 1000) == [(0, 40), (299, 330), (499, 580), (959, 1000)]


def test_vector_regions_far_apart_stay_separate():
    assert vecscreen._vector_regions([(400, 420), (100, 150)], 1000) == [(99, 150), (399, 420)]


# --- _write_cleaned -----------------------------------------------------------------------------


def test_terminal_trim_and_internal_split(tmp_path):
    fa = tmp_path / "in.fa"
    seq_t = "N" * 20 + "A" * 300 + "T" * 30  # 350 bp: vector 11..20 (extended to 1) and 321..350
    seq_s = "A" * 250 + "G" * 50 + "C" * 250  # internal vector 251..300
    fa.write_text(f">term\n{seq_t}\n>split\n{seq_s}\n>short\n{'A' * 199}\n>ok\n{'A' * 200}\n")
    hits = {"term": [(11, 20), (321, 350)], "split": [(251, 300), (251, 300)]}
    out = tmp_path / "out.fa"
    assert vecscreen._write_cleaned(str(fa), hits, str(out)) == (4, 4)
    assert _read_fasta(out) == {"term": "A" * 300, "split1_split": "A" * 250, "split2_split": "C" * 250, "ok": "A" * 200}


def test_pieces_shorter_than_200_dropped(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">c\n" + "A" * 100 + "G" * 50 + "C" * 250 + "\n>t\n" + "A" * 230 + "\n")
    out = tmp_path / "out.fa"
    vecscreen._write_cleaned(str(fa), {"c": [(101, 150)], "t": [(201, 230)]}, str(out))
    assert _read_fasta(out) == {"split2_c": "C" * 250, "t": "A" * 200}  # 100 bp piece dropped; 200 bp kept


def test_mostly_vector_contig_excluded(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">v\n" + "A" * 1000 + "\n>half\n" + "A" * 1000 + "\n")
    out = tmp_path / "out.fa"
    vecscreen._write_cleaned(str(fa), {"v": [(100, 300), (340, 700)], "half": [(1, 500)]}, str(out))
    assert _read_fasta(out) == {"half": "A" * 500}  # exactly 50% is kept, > 50% excluded


# --- run() --------------------------------------------------------------------------------------


@pytest.mark.parametrize("pipe", [False, True])
def test_run_end_to_end(tmp_path, monkeypatch, caplog, pipe):
    monkeypatch.chdir(tmp_path)
    infile = tmp_path / "asm.fasta"
    infile.write_text(">keep\n" + "A" * 300 + "\n>trim\n" + "C" * 300 + "\n")
    cmds = []

    def fake_run(cmd, *args, **kw):
        cmds.append(cmd)
        Path(cmd[cmd.index("-out") + 1]).write_text(_vec_row("trim", 1, 30, 25, qlen=300) + "\n" + _vec_row("keep", 120, 150, 29, qlen=300) + "\n")
        return Namespace(returncode=0)

    caplog.set_level(logging.INFO)
    (tmp_path / "sub").mkdir()
    with patch("aaftf.vecscreen.require_tools"), patch("aaftf.vecscreen.require_contamdb", return_value="/db/contamdb"), patch("aaftf.vecscreen.run_cmd", side_effect=fake_run):
        vecscreen.run(infile=str(infile), outfile="sub/asm.vecscreen.fasta", workdir=str(tmp_path / "wd"), pipe=pipe)
    (cmd,) = cmds
    assert cmd[: 1 + len(vecscreen.VECSCREEN_TASK_ARGS)] == ["blastn", *vecscreen.VECSCREEN_TASK_ARGS] and cmd[cmd.index("-db") + 1] == "/db/contamdb"
    assert _read_fasta(tmp_path / "sub" / "asm.vecscreen.fasta") == {"keep": "A" * 300, "trim": "C" * 270}
    assert ("AAFTF taxpurge" in caplog.text) is (not pipe)


def test_run_needs_the_database(tmp_path, monkeypatch):
    monkeypatch.setenv("AAFTF_DB", str(tmp_path))
    with patch("aaftf.vecscreen.require_tools"), pytest.raises(FileNotFoundError, match="AAFTF database contamdb"):
        vecscreen.run(infile="in.fa", outfile="out.fa")
