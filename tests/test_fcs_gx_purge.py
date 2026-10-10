"""Unit tests for aaftf/fcs_gx_purge.py (run_gx.py is mocked)."""

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.fcs_gx_purge import run

pytestmark = pytest.mark.unit

_HEADER = "##[header]\n#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n"


def _run(tmp_path, report_rows, returncode=0, write_report=True):
    """Run fcs_gx_purge.run() on a 3-contig assembly with a fake run_gx.py; return the kept contig IDs."""
    asm = tmp_path / "asm.fasta"
    asm.write_text(">keep\nACGT\n>bad\nACGT\n>trimme\nACGT\n")
    db = tmp_path / "gxdb" / "all"
    db.parent.mkdir()
    Path(f"{db}.gxi").write_text("")
    workdir = tmp_path / "wd"
    out = tmp_path / "out.fasta"

    def _fake_run(cmd, **kwargs):
        if write_report:
            Path(workdir, "asm.4890.fcs_gx_report.txt").write_text(_HEADER + "".join("\t".join(r) + "\n" for r in report_rows))
        return subprocess.CompletedProcess(cmd, returncode)

    with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run):
        run(input=str(asm), outfile=str(out), workdir=str(workdir), db=str(db))
    return [line[1:] for line in out.read_text().split() if line.startswith(">")]


def test_only_exclude_rows_are_dropped(tmp_path, caplog):
    rows = [["bad", "1", "4", "4", "EXCLUDE", "anml:insects", "100", "fly"], ["trimme", "1", "2", "4", "TRIM", "prok:bact", "50", "E. coli"]]
    assert _run(tmp_path, rows) == ["keep", "trimme"]
    assert "TRIM for trimme:1..2" in caplog.text


def test_action_prefix_accepted(tmp_path):
    assert _run(tmp_path, [["bad", "1", "4", "4", "ACTION_EXCLUDE", "x", "100", "y"]]) == ["keep", "trimme"]


def test_failed_run_gx_raises(tmp_path):
    with pytest.raises(RuntimeError, match="run_gx.py failed"):
        _run(tmp_path, [], returncode=1)


def test_missing_report_raises(tmp_path):
    with pytest.raises(RuntimeError, match="did not write"):
        _run(tmp_path, [], write_report=False)


class TestNextCommandHint:
    """The "next command" hint is logged unless -q/--quiet is given."""

    def _main(self, tmp_path, *extra):
        import sys

        from aaftf.main import main

        asm = tmp_path / "asm.fasta"
        asm.write_text(">keep\nACGT\n")
        db = tmp_path / "gxdb" / "all"
        db.parent.mkdir()
        Path(f"{db}.gxi").write_text("")
        workdir = tmp_path / "wd"

        def _fake_run(cmd, **kwargs):
            Path(workdir, "asm.4890.fcs_gx_report.txt").write_text(_HEADER)
            return subprocess.CompletedProcess(cmd, 0)

        argv = ["AAFTF", "fcs_gx_purge", "-i", str(asm), "-o", str(tmp_path / "out.fasta"), "-w", str(workdir), "-d", str(db), *extra]
        with patch.object(sys, "argv", argv), patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run):
            assert main() == 0

    def test_hint_shown_by_default(self, tmp_path, capsys):
        self._main(tmp_path)
        assert "Your next command might be:" in capsys.readouterr().err

    def test_hint_hidden_with_quiet(self, tmp_path, capsys):
        self._main(tmp_path, "-q")
        assert "Your next command might be:" not in capsys.readouterr().err


@pytest.mark.parametrize("pipe, shown", [(False, True), (True, False)])
def test_pipe_controls_hint(tmp_path, caplog, pipe, shown):
    """run(pipe=True), as pipeline calls it, logs no "next command" hint (pipe is not a CLI option)."""
    import logging

    asm = tmp_path / "asm.fasta"
    asm.write_text(">keep\nACGT\n")
    db = tmp_path / "gxdb" / "all"
    db.parent.mkdir()
    Path(f"{db}.gxi").write_text("")
    workdir = tmp_path / "wd"

    def _fake_run(cmd, **kwargs):
        Path(workdir, "asm.4890.fcs_gx_report.txt").write_text(_HEADER)
        return subprocess.CompletedProcess(cmd, 0)

    with caplog.at_level(logging.INFO, logger="aaftf"), patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run):
        run(input=str(asm), outfile=str(tmp_path / "out.fasta"), workdir=str(workdir), db=str(db), pipe=pipe)
    assert ("Your next command might be:" in caplog.text) is shown
