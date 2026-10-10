"""Per-step log files: <workdir>/<step>.log, kept exactly when the work directory is kept."""

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.main import main

pytestmark = pytest.mark.unit

_REPORT_HEADER = "#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n"


def _fcs_gx(tmp_path, *extra, run_gx_rc=0):
    """Run `AAFTF fcs_gx_purge` through main() with run_gx.py mocked; return main()'s exit code."""
    asm = tmp_path / "asm.fasta"
    asm.write_text(">keep\nACGT\n")
    db = tmp_path / "gxdb" / "all"
    db.parent.mkdir(exist_ok=True)
    Path(f"{db}.gxi").write_text("")

    def _fake_run(cmd, **kwargs):
        outdir = cmd[cmd.index("--out-dir") + 1]
        Path(outdir, "asm.4890.fcs_gx_report.txt").write_text(_REPORT_HEADER)
        return subprocess.CompletedProcess(cmd, run_gx_rc)

    argv = ["AAFTF", "fcs_gx_purge", "-i", str(asm), "-o", str(tmp_path / "out.fasta"), "-d", str(db), *extra]
    with patch.object(sys, "argv", argv), patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run):
        return main()


class TestWorkdirLog:
    def test_own_workdir_keeps_full_log(self, tmp_path):
        assert _fcs_gx(tmp_path, "-w", str(tmp_path / "wd")) == 0
        log = (tmp_path / "wd" / "fcs_gx_purge.log").read_text()
        assert "Running AAFTF" in log  # logged before the work directory existed
        assert "fcs-gx assembly is" in log and "Your next command might be:" in log
        assert "\033[" not in log  # no terminal colour codes

    def test_quiet_still_logs_info_to_file(self, tmp_path, capsys):
        assert _fcs_gx(tmp_path, "-w", str(tmp_path / "wd"), "-q") == 0
        assert "Running AAFTF" not in capsys.readouterr().err
        assert "Running AAFTF" in (tmp_path / "wd" / "fcs_gx_purge.log").read_text()

    def test_reruns_append(self, tmp_path):
        _fcs_gx(tmp_path, "-w", str(tmp_path / "wd"))
        _fcs_gx(tmp_path, "-w", str(tmp_path / "wd"))
        assert (tmp_path / "wd" / "fcs_gx_purge.log").read_text().count("Running AAFTF") == 2

    def test_auto_workdir_and_its_log_removed(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert _fcs_gx(tmp_path) == 0
        assert not list(tmp_path.glob("aaftf-fcs_gx_purge_*")) and not list(tmp_path.glob("*.log"))

    def test_auto_workdir_log_kept_with_verbose(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert _fcs_gx(tmp_path, "-v") == 0
        (workdir,) = tmp_path.glob("aaftf-fcs_gx_purge_*")
        assert "Running AAFTF" in (workdir / "fcs_gx_purge.log").read_text()

    def test_failure_keeps_log_with_traceback(self, tmp_path, capsys):
        assert _fcs_gx(tmp_path, "-w", str(tmp_path / "wd"), run_gx_rc=1) == 1
        log = (tmp_path / "wd" / "fcs_gx_purge.log").read_text()
        assert "ERROR: An error occurred: run_gx.py failed" in log and "Traceback (most recent call last)" in log
        assert "Traceback" not in capsys.readouterr().err  # the terminal shows it only with -v


class TestNoWorkdirLog:
    """Subcommands without a work directory (e.g. sort) write ./<command>.log only with -v."""

    def _sort(self, tmp_path, monkeypatch, *extra):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "in.fa").write_text(">a\nACGT\n")
        with patch.object(sys, "argv", ["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", *extra]):
            assert main() == 0

    def test_verbose_writes_command_log(self, tmp_path, monkeypatch):
        self._sort(tmp_path, monkeypatch, "-v")
        assert "Running AAFTF" in (tmp_path / "sort.log").read_text()

    def test_no_log_without_verbose(self, tmp_path, monkeypatch):
        self._sort(tmp_path, monkeypatch)
        assert not (tmp_path / "sort.log").exists()


def test_two_steps_like_pipeline(tmp_path, monkeypatch):
    """Each step log holds its own messages; between-step lines go to the next step, the final hint to the last."""
    import logging

    from aaftf import utility

    monkeypatch.chdir(tmp_path)
    log = logging.getLogger("aaftf.pipeline")
    utility.setup_logging()
    dir_a, custom_a = utility.make_workdir(str(tmp_path / "a"), "stepa")
    log.info("inside stepa")
    utility.cleanup_workdir(dir_a, False, custom_a)
    log.info("between steps")
    dir_b, custom_b = utility.make_workdir(str(tmp_path / "b"), "stepb")
    log.info("inside stepb")
    utility.cleanup_workdir(dir_b, False, custom_b)
    log.info("final hint")
    utility.finish_logging("pipeline", False)

    log_a = (tmp_path / "a" / "stepa.log").read_text()
    log_b = (tmp_path / "b" / "stepb.log").read_text()
    assert "inside stepa" in log_a
    assert not any(msg in log_a for msg in ("between steps", "inside stepb", "final hint"))
    assert "inside stepa" not in log_b
    assert log_b.index("between steps") < log_b.index("inside stepb") < log_b.index("final hint")
    assert log_b.count("final hint") == 1
    assert not (tmp_path / "pipeline.log").exists()
