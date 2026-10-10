"""Unit tests for aaftf/polish.py.

Covers:
  - CLI parser defaults and flag presence
  - run() guard conditions (no reads, racon/nextpolish2 without longreads)
  - pypolca failure/success (files copied to expected destinations)
  - polypolish failure/success
  - nextpolish2 failure/success

No external bioinformatics tools are invoked in the unit tests.
"""

import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import aaftf.polish as polish
from aaftf.main import main

# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF polish ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_polish(argv):
    captured = {}

    def _capture(**kwargs):
        captured["args"] = Namespace(**kwargs)

    with patch.object(sys, "argv", argv):
        with patch("aaftf.polish.run", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for polish.run()
# ---------------------------------------------------------------------------


def _make_args(tmp_path, method="pypolca", **overrides):
    defaults = dict(
        method=method,
        memory=4,
        cpus=2,
        read1=str(tmp_path / "R1.fq"),
        read2=str(tmp_path / "R2.fq"),
        longreads=None,
        workdir=str(tmp_path / "workdir"),
        infile=str(tmp_path / "asm.fa"),
        outfile=None,
        debug=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


# ---------------------------------------------------------------------------
# Parser defaults and required flags
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestPolishParser:
    def test_debug_false_by_default(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq", "-v"])
        assert args.debug is True

    def test_default_method_is_polypolish(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq"])
        assert args.method == "polypolish"

    def test_method_polypolish(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq", "--method", "polypolish"])
        assert args.method == "polypolish"

    def test_method_nextpolish2(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-lr", "lr.fq", "--method", "nextpolish2"])
        assert args.method == "nextpolish2"

    def test_method_racon(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-lr", "lr.fq", "--method", "racon"])
        assert args.method == "racon"

    def test_invalid_method_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "polish", "-i", "asm.fa", "--method", "pilon"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_default_cpus(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq", "-c", "8"])
        assert args.cpus == 8

    def test_default_memory(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq"])
        assert args.memory == 16

    def test_custom_memory(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq", "-m", "32"])
        assert args.memory == 32

    def test_parses_read1_reads(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq"])
        assert args.read1 == "R1.fq"

    def test_parses_read2_reads(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-1", "R1.fq", "-2", "R2.fq"])
        assert args.read2 == "R2.fq"

    def test_parses_longreads(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-lr", "lr.fq"])
        assert args.longreads == "lr.fq"

    def test_missing_infile_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "polish", "-1", "R1.fq"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_polish_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "polish", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0


# ---------------------------------------------------------------------------
# run() guard conditions — exit before any external tool is called
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestPolishRunGuards:
    def test_racon_without_longreads_raises(self, tmp_path):
        args = _make_args(tmp_path, method="racon", longreads=None, read1=None, read2=None)
        from aaftf.polish import run

        with pytest.raises(ValueError):
            run(**vars(args))

    def test_nextpolish2_without_longreads_raises(self, tmp_path):
        args = _make_args(tmp_path, method="nextpolish2", longreads=None)
        from aaftf.polish import run

        with pytest.raises(ValueError):
            run(**vars(args))

    def test_pypolca_without_reads_raises(self, tmp_path):
        args = _make_args(tmp_path, method="pypolca", read1=None, read2=None)
        from aaftf.polish import run

        with pytest.raises(ValueError):
            run(**vars(args))

    def test_polypolish_without_reads_raises(self, tmp_path):
        args = _make_args(tmp_path, method="polypolish", read1=None, read2=None)
        from aaftf.polish import run

        with pytest.raises(ValueError):
            run(**vars(args))

    def test_polypolish_without_read2_reads_raises(self, tmp_path):
        """Polypolish requires paired reads."""
        (tmp_path / "R1.fq").write_text("@r\nA\n+\nI\n")
        args = _make_args(tmp_path, method="polypolish", read2=None)
        from aaftf.polish import run

        with patch("aaftf.polish.shutil.which", return_value=True):
            with pytest.raises(ValueError):
                run(**vars(args))


# ---------------------------------------------------------------------------
# pypolca failure/success detection
# ---------------------------------------------------------------------------


def _setup_pypolca_run(tmp_path, outfile=None):
    """Create input files and return a pypolca-method Namespace."""
    infile = tmp_path / "asm.fa"
    infile.write_text(">contig1\nATCGATCG\n")
    workdir = tmp_path / "wdir"
    workdir.mkdir()
    r1 = tmp_path / "R1.fq"
    r2 = tmp_path / "R2.fq"
    r1.write_text("@r\nA\n+\nI\n")
    r2.write_text("@r\nA\n+\nI\n")
    return _make_args(
        tmp_path,
        method="pypolca",
        infile=str(infile),
        workdir=str(workdir),
        read1=str(r1),
        read2=str(r2),
        outfile=outfile,
    )


@pytest.mark.unit
class TestPolishPypolca:
    def test_nonzero_exit_raises_runtime_error(self, tmp_path):
        args = _setup_pypolca_run(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 1

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", return_value=mock_result):
            with patch("aaftf.polish.shutil.which", return_value=True):
                with pytest.raises(RuntimeError):
                    run(**vars(args))

    def test_missing_output_file_raises_runtime_error(self, tmp_path):
        """pypolca exits 0 but creates no output file — must still raise."""
        args = _setup_pypolca_run(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 0
        # Output file deliberately NOT created

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", return_value=mock_result):
            with patch("aaftf.polish.shutil.which", return_value=True):
                with pytest.raises(RuntimeError):
                    run(**vars(args))

    def test_success_copies_corrected_fasta(self, tmp_path):
        """pypolca success: {prefix}_corrected.fasta is copied to the outfile path."""
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_pypolca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0

        def _fake_run(cmd, **kw):
            outdir = workdir / "pypolca_out"
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "pypolca_corrected.fasta").write_text(">contig1\nATCGATCG\n")
            (outdir / "pypolca.vcf").write_text("")
            (outdir / "pypolca.report").write_text("")
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        assert (tmp_path / "polished.fasta").exists()

    def test_success_copies_vcf(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_pypolca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0

        def _fake_run(cmd, **kw):
            outdir = workdir / "pypolca_out"
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "pypolca_corrected.fasta").write_text(">contig1\nATCGATCG\n")
            (outdir / "pypolca.vcf").write_text("##vcf\n")
            (outdir / "pypolca.report").write_text("")
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        assert (tmp_path / "polished.fasta.vcf").exists()

    def test_success_copies_report(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_pypolca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0

        def _fake_run(cmd, **kw):
            outdir = workdir / "pypolca_out"
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "pypolca_corrected.fasta").write_text(">contig1\nATCGATCG\n")
            (outdir / "pypolca.vcf").write_text("")
            (outdir / "pypolca.report").write_text("POLCA report\n")
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        assert (tmp_path / "polished.fasta.pypolca_report.txt").exists()

    def test_pypolca_cmd_includes_reads(self, tmp_path):
        """pypolca subprocess call must include the read files."""
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_pypolca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0
        captured_cmds = []

        def _fake_run(cmd, **kw):
            captured_cmds.append(cmd)
            outdir = workdir / "pypolca_out"
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "pypolca_corrected.fasta").write_text(">c\nATCG\n")
            (outdir / "pypolca.vcf").write_text("")
            (outdir / "pypolca.report").write_text("")
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        pypolca_cmd = captured_cmds[0]
        r1 = next((pypolca_cmd[i + 1] for i, a in enumerate(pypolca_cmd) if a == "-1"), None)
        r2 = next((pypolca_cmd[i + 1] for i, a in enumerate(pypolca_cmd) if a == "-2"), None)
        assert r1 is not None and "R1.fq" in r1
        assert r2 is not None and "R2.fq" in r2


# ---------------------------------------------------------------------------
# polypolish failure/success detection
# ---------------------------------------------------------------------------


def _setup_polypolish_run(tmp_path, outfile=None):
    """Create input files and return a polypolish-method Namespace."""
    infile = tmp_path / "asm.fa"
    infile.write_text(">contig1\nATCGATCG\n")
    workdir = tmp_path / "wdir"
    workdir.mkdir()
    r1 = tmp_path / "R1.fq"
    r2 = tmp_path / "R2.fq"
    r1.write_text("@r\nA\n+\nI\n")
    r2.write_text("@r\nA\n+\nI\n")
    return _make_args(
        tmp_path,
        method="polypolish",
        infile=str(infile),
        workdir=str(workdir),
        read1=str(r1),
        read2=str(r2),
        outfile=outfile,
    )


@pytest.mark.unit
class TestPolishPolypolish:
    def test_empty_output_raises_runtime_error(self, tmp_path):
        """polypolish polish writes to stdout; an empty result must still raise."""
        args = _setup_polypolish_run(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 0

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", return_value=mock_result):
            with patch("aaftf.polish.shutil.which", return_value=True):
                with pytest.raises(RuntimeError):
                    run(**vars(args))

    def test_success_copies_corrected_fasta(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_polypolish_run(tmp_path, outfile=outfile)

        def _fake_run(cmd, stdout=None, **kw):
            if cmd[:2] == ["polypolish", "polish"] and stdout is not None:
                stdout.write(">contig1\nATCGATCG\n")
            mock_result = MagicMock()
            mock_result.returncode = 0
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        assert (tmp_path / "polished.fasta").exists()
        assert (tmp_path / "polished.fasta").stat().st_size > 0

    def test_bwa_mem_cmd_includes_reads(self, tmp_path):
        """Each read file is aligned separately with `bwa mem -a`."""
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_polypolish_run(tmp_path, outfile=outfile)
        captured_cmds = []

        def _fake_run(cmd, stdout=None, **kw):
            captured_cmds.append(cmd)
            if cmd[:2] == ["polypolish", "polish"] and stdout is not None:
                stdout.write(">contig1\nATCGATCG\n")
            mock_result = MagicMock()
            mock_result.returncode = 0
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        bwa_mem_cmds = [c for c in captured_cmds if c[:2] == ["bwa", "mem"]]
        assert len(bwa_mem_cmds) == 2
        assert all("-a" in c for c in bwa_mem_cmds)
        reads_used = {c[-1] for c in bwa_mem_cmds}
        assert any("R1.fq" in r for r in reads_used)
        assert any("R2.fq" in r for r in reads_used)


# ---------------------------------------------------------------------------
# racon failure/success detection
# ---------------------------------------------------------------------------


def _setup_racon_run(tmp_path, outfile=None):
    """Create input files and return a racon-method Namespace."""
    infile = tmp_path / "asm.fa"
    infile.write_text(">contig1\nATCGATCG\n")
    workdir = tmp_path / "wdir"
    workdir.mkdir()
    lr = tmp_path / "lr.fq"
    lr.write_text("@r\nA\n+\nI\n")
    return _make_args(
        tmp_path,
        method="racon",
        infile=str(infile),
        workdir=str(workdir),
        read1=None,
        read2=None,
        longreads=str(lr),
        outfile=outfile,
    )


@pytest.mark.unit
class TestPolishRacon:
    def test_empty_output_raises_runtime_error(self, tmp_path):
        """racon writes to stdout; an empty result must still raise."""
        args = _setup_racon_run(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 0

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", return_value=mock_result):
            with patch("aaftf.polish.shutil.which", return_value=True):
                with pytest.raises(RuntimeError):
                    run(**vars(args))

    def test_success_copies_corrected_fasta(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_racon_run(tmp_path, outfile=outfile)

        def _fake_run(cmd, stdout=None, **kw):
            if cmd[:1] == ["racon"] and stdout is not None:
                stdout.write(">contig1\nATCGATCG\n")
            mock_result = MagicMock()
            mock_result.returncode = 0
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        assert (tmp_path / "polished.fasta").exists()
        assert (tmp_path / "polished.fasta").stat().st_size > 0

    def test_racon_cmd_includes_longreads_and_overlaps(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_racon_run(tmp_path, outfile=outfile)
        captured_cmds = []

        def _fake_run(cmd, stdout=None, **kw):
            captured_cmds.append(cmd)
            if cmd[:1] == ["racon"] and stdout is not None:
                stdout.write(">contig1\nATCGATCG\n")
            mock_result = MagicMock()
            mock_result.returncode = 0
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.shutil.which", return_value=True):
                run(**vars(args))

        minimap_cmds = [c for c in captured_cmds if c[:1] == ["minimap2"]]
        racon_cmds = [c for c in captured_cmds if c[:1] == ["racon"]]
        assert len(minimap_cmds) == 1
        assert len(racon_cmds) == 1
        assert any("lr.fq" in a for a in minimap_cmds[0])
        assert any("lr.fq" in a for a in racon_cmds[0])
        assert "overlaps.paf" in racon_cmds[0]


# ---------------------------------------------------------------------------
# nextpolish2 failure/success detection
# ---------------------------------------------------------------------------


def _setup_nextpolish2_run(tmp_path, outfile=None):
    """Create input files and return a nextpolish2-method Namespace."""
    infile = tmp_path / "asm.fa"
    infile.write_text(">contig1\nATCGATCG\n")
    workdir = tmp_path / "wdir"
    workdir.mkdir()
    r1 = tmp_path / "R1.fq"
    r2 = tmp_path / "R2.fq"
    lr = tmp_path / "lr.fq"
    r1.write_text("@r\nA\n+\nI\n")
    r2.write_text("@r\nA\n+\nI\n")
    lr.write_text("@r\nA\n+\nI\n")
    return _make_args(
        tmp_path,
        method="nextpolish2",
        infile=str(infile),
        workdir=str(workdir),
        read1=str(r1),
        read2=str(r2),
        longreads=str(lr),
        outfile=outfile,
    )


@pytest.mark.unit
class TestPolishNextpolish2:
    def test_nonzero_exit_raises_runtime_error(self, tmp_path):
        args = _setup_nextpolish2_run(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_popen = MagicMock()
        mock_popen.stdout = MagicMock()
        mock_popen.wait.return_value = 0
        mock_popen.returncode = 0

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", return_value=mock_result):
            with patch("aaftf.polish.subprocess.Popen", return_value=mock_popen):
                with patch("aaftf.polish.shutil.which", return_value=True):
                    with pytest.raises(RuntimeError):
                        run(**vars(args))

    def test_success_copies_corrected_fasta(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_nextpolish2_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_popen = MagicMock()
        mock_popen.stdout = MagicMock()
        mock_popen.wait.return_value = 0
        mock_popen.returncode = 0

        def _fake_run(cmd, **kw):
            if cmd[0] == "nextPolish2":
                out_idx = cmd.index("-o") + 1
                (workdir / cmd[out_idx]).write_text(">contig1\nATCGATCG\n")
            mock_result = MagicMock()
            mock_result.returncode = 0
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.subprocess.Popen", return_value=mock_popen):
                with patch("aaftf.polish.shutil.which", return_value=True):
                    run(**vars(args))

        assert (tmp_path / "polished.fasta").exists()

    def test_yak_count_includes_short_reads(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_nextpolish2_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_popen = MagicMock()
        mock_popen.stdout = MagicMock()
        mock_popen.wait.return_value = 0
        mock_popen.returncode = 0
        captured_cmds = []

        def _fake_run(cmd, **kw):
            captured_cmds.append(cmd)
            if cmd[0] == "nextPolish2":
                out_idx = cmd.index("-o") + 1
                (workdir / cmd[out_idx]).write_text(">contig1\nATCGATCG\n")
            mock_result = MagicMock()
            mock_result.returncode = 0
            return mock_result

        from aaftf.polish import run

        with patch("aaftf.polish.subprocess.run", side_effect=_fake_run):
            with patch("aaftf.polish.subprocess.Popen", return_value=mock_popen):
                with patch("aaftf.polish.shutil.which", return_value=True):
                    run(**vars(args))

        yak_cmds = [c for c in captured_cmds if c[:2] == ["yak", "count"]]
        assert len(yak_cmds) == 1
        assert any("R1.fq" in a for a in yak_cmds[0])
        assert any("R2.fq" in a for a in yak_cmds[0])


# ---------------------------------------------------------------------------
# polca (aaftf_polca.sh) threads/memory/environment
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestPolishPolca:
    def test_threads_and_memory_passed_to_polca(self, tmp_path, monkeypatch):
        """cpus/memory reach aaftf_polca.sh as -t/-m and NUM_THREADS/MEM; BASH_FUNC_* is dropped."""
        (tmp_path / "asm.fa").write_text(">contig1\nATCG\n")
        workdir = tmp_path / "workdir"
        workdir.mkdir()
        monkeypatch.setenv("BASH_FUNC_which%%", "() { :; }")
        calls = {}

        def fake_run(cmd, cwd=None, env=None, **kw):
            calls.update(cmd=cmd, env=env)
            Path(cwd, "asm.fa.PolcaCorrected.fa").write_text(">contig1\nATCG\n")
            Path(cwd, "asm.fa.report").write_text("report\n")
            return MagicMock(returncode=0)

        with patch("aaftf.polish.subprocess.run", side_effect=fake_run):
            ret, out = polish.run_polca(str(tmp_path / "asm.fa"), "R1.fq", "R2.fq", 8, 32, str(workdir), "polca.log", str(tmp_path / "out.fasta"))
        assert ret == 0 and out.endswith("asm.fa.PolcaCorrected.fa")
        cmd = calls["cmd"]
        assert cmd[cmd.index("-t") + 1] == "8" and cmd[cmd.index("-m") + 1] == "4G"
        assert cmd[cmd.index("-r") + 1] == "R1.fq R2.fq"
        assert calls["env"]["NUM_THREADS"] == "8" and calls["env"]["MEM"] == "4G"
        assert not any(k.startswith("BASH_FUNC_") for k in calls["env"])
        assert (tmp_path / "out.fasta.polca_report.txt").exists()
