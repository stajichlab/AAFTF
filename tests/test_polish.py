"""Unit tests for AAFTF/polish.py.

Covers:
  - CLI parser defaults and flag presence
  - run() guard conditions (no reads, zero iterations, racon without longreads)
  - polca failure detection (non-zero exit code, missing output file)
  - polca success path (files copied to expected destinations)
  - pilon convergence (stops at iteration where zero changes are found)
  - Integration test: real pilon run on Rhizopus test data (requires pilon + bwa)

No external bioinformatics tools are invoked in the unit tests.
"""

import shutil
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from AAFTF.AAFTF_main import main

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF polish ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_polish(argv):
    captured = {}

    def _capture(parser, args):
        captured["args"] = args

    with patch.object(sys, "argv", argv):
        with patch("AAFTF.AAFTF_main.run_subtool", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for polish.run()
# ---------------------------------------------------------------------------


def _make_args(tmp_path, method="pilon", **overrides):
    defaults = dict(
        method=method,
        memory=4,
        cpus=2,
        left=str(tmp_path / "R1.fq"),
        right=str(tmp_path / "R2.fq"),
        longreads=None,
        workdir=str(tmp_path / "workdir"),
        infile=str(tmp_path / "asm.fa"),
        outfile=None,
        iterations=1,
        debug=False,
        diploid=False,
        ploidy=1,
        pipe=True,
        polca="polca.sh",
        polca_samtools=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


# ---------------------------------------------------------------------------
# Parser defaults and required flags
# ---------------------------------------------------------------------------


class TestPolishParser:
    def test_debug_false_by_default(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "-v"])
        assert args.debug is True

    def test_pipe_false_by_default(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.pipe is False

    def test_pipe_flag_sets_true(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "--pipe"])
        assert args.pipe is True

    def test_default_method_is_pilon(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.method == "pilon"

    def test_method_polca(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "--method", "polca"])
        assert args.method == "polca"

    def test_method_nextpolish(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "--method", "nextpolish"])
        assert args.method == "nextpolish"

    def test_method_racon(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-lr", "lr.fq", "--method", "racon"])
        assert args.method == "racon"

    def test_default_iterations(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.iterations == 5

    def test_custom_iterations(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "-it", "2"])
        assert args.iterations == 2

    def test_default_cpus(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "-c", "8"])
        assert args.cpus == 8

    def test_default_memory(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.memory == 16

    def test_custom_memory(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "-m", "32"])
        assert args.memory == 32

    def test_default_polca_exe(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.polca == "polca.sh"

    def test_custom_polca_exe(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "--polca", "/opt/bin/polca.sh"])
        assert args.polca == "/opt/bin/polca.sh"

    def test_diploid_false_by_default(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.diploid is False

    def test_diploid_flag_sets_true(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "--diploid"])
        assert args.diploid is True

    def test_default_ploidy(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.ploidy == 1

    def test_parses_left_reads(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.left == "R1.fq"

    def test_parses_right_reads(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "-r", "R2.fq"])
        assert args.right == "R2.fq"

    def test_parses_longreads(self):
        args = _parse_polish(["AAFTF", "polish", "-i", "asm.fa", "-lr", "lr.fq"])
        assert args.longreads == "lr.fq"

    def test_pilon_alias_routes_to_polish(self):
        args = _parse_polish(["AAFTF", "pilon", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.infile == "asm.fa"

    def test_polca_alias_routes_to_polish(self):
        args = _parse_polish(["AAFTF", "polca", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.infile == "asm.fa"

    def test_missing_infile_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "polish", "-l", "R1.fq"]):
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


class TestPolishRunGuards:
    def test_racon_without_longreads_exits(self, tmp_path):
        args = _make_args(tmp_path, method="racon", longreads=None, left=None, right=None)
        from AAFTF.polish import run

        with pytest.raises(SystemExit):
            run(None, args)

    def test_pilon_without_reads_exits(self, tmp_path):
        args = _make_args(tmp_path, method="pilon", left=None, right=None)
        from AAFTF.polish import run

        with pytest.raises(SystemExit):
            run(None, args)

    def test_polca_without_reads_exits(self, tmp_path):
        args = _make_args(tmp_path, method="polca", left=None, right=None)
        from AAFTF.polish import run

        with pytest.raises(SystemExit):
            run(None, args)

    def test_zero_iterations_exits(self, tmp_path):
        (tmp_path / "R1.fq").write_text("@r\nA\n+\nI\n")
        (tmp_path / "R2.fq").write_text("@r\nA\n+\nI\n")
        args = _make_args(tmp_path, method="pilon", iterations=0)
        from AAFTF.polish import run

        with pytest.raises(SystemExit):
            run(None, args)

    def test_negative_iterations_exits(self, tmp_path):
        (tmp_path / "R1.fq").write_text("@r\nA\n+\nI\n")
        (tmp_path / "R2.fq").write_text("@r\nA\n+\nI\n")
        args = _make_args(tmp_path, method="pilon", iterations=-2)
        from AAFTF.polish import run

        with pytest.raises(SystemExit):
            run(None, args)


# ---------------------------------------------------------------------------
# polca failure detection
# ---------------------------------------------------------------------------


def _setup_polca_run(tmp_path, outfile=None):
    """Create input files and return a polca-method Namespace."""
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
        method="polca",
        infile=str(infile),
        workdir=str(workdir),
        left=str(r1),
        right=str(r2),
        outfile=outfile,
    )


class TestPolishPolcaFailure:
    def test_nonzero_exit_raises_systemexit(self, tmp_path):
        args = _setup_polca_run(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 1

        from packaging.version import Version

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", return_value=mock_result):
            with patch("AAFTF.polish.get_samtools_version", return_value=Version("1.23")):
                with pytest.raises(SystemExit) as exc:
                    run(None, args)
        assert exc.value.code == 1

    def test_missing_output_file_raises_systemexit(self, tmp_path):
        """polca exits 0 but creates no output file — must still sys.exit(1)."""
        args = _setup_polca_run(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 0
        # Output file deliberately NOT created

        from packaging.version import Version

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", return_value=mock_result):
            with patch("AAFTF.polish.get_samtools_version", return_value=Version("1.23")):
                with pytest.raises(SystemExit) as exc:
                    run(None, args)
        assert exc.value.code == 1

    def test_success_copies_corrected_fasta(self, tmp_path):
        """polca success: PolcaCorrected.fa is copied to the outfile path."""
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_polca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0

        def _fake_run(cmd, **kw):
            (workdir / "asm.fa.PolcaCorrected.fa").write_text(">contig1\nATCGATCG\n")
            (workdir / "asm.fa.vcf").write_text("")
            (workdir / "asm.fa.report").write_text("")
            return mock_result

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", side_effect=_fake_run):
            run(None, args)

        assert (tmp_path / "polished.fasta").exists()

    def test_success_copies_vcf(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_polca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0

        def _fake_run(cmd, **kw):
            (workdir / "asm.fa.PolcaCorrected.fa").write_text(">contig1\nATCGATCG\n")
            (workdir / "asm.fa.vcf").write_text("##vcf\n")
            (workdir / "asm.fa.report").write_text("")
            return mock_result

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", side_effect=_fake_run):
            run(None, args)

        assert (tmp_path / "polished.fasta.vcf").exists()

    def test_success_copies_report(self, tmp_path):
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_polca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0

        def _fake_run(cmd, **kw):
            (workdir / "asm.fa.PolcaCorrected.fa").write_text(">contig1\nATCGATCG\n")
            (workdir / "asm.fa.vcf").write_text("")
            (workdir / "asm.fa.report").write_text("POLCA report\n")
            return mock_result

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", side_effect=_fake_run):
            run(None, args)

        assert (tmp_path / "polished.fasta.polca_report.txt").exists()

    def test_polca_cmd_includes_reads(self, tmp_path):
        """polca subprocess call must include the read files."""
        outfile = str(tmp_path / "polished.fasta")
        args = _setup_polca_run(tmp_path, outfile=outfile)
        workdir = tmp_path / "wdir"
        mock_result = MagicMock()
        mock_result.returncode = 0
        captured_cmds = []

        def _fake_run(cmd, **kw):
            captured_cmds.append(cmd)
            (workdir / "asm.fa.PolcaCorrected.fa").write_text(">c\nATCG\n")
            (workdir / "asm.fa.vcf").write_text("")
            (workdir / "asm.fa.report").write_text("")
            return mock_result

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", side_effect=_fake_run):
            run(None, args)

        polca_cmd = captured_cmds[0]
        reads_arg = next((polca_cmd[i + 1] for i, a in enumerate(polca_cmd) if a == "-r"), None)
        assert reads_arg is not None
        assert "R1.fq" in reads_arg
        assert "R2.fq" in reads_arg


# ---------------------------------------------------------------------------
# pilon: convergence (zero changes stops the loop)
# ---------------------------------------------------------------------------


class TestPolishPilonConvergence:
    def test_stops_after_zero_changes(self, tmp_path):
        """Pilon loop breaks at iteration 1 when the .changes file is empty."""
        infile = tmp_path / "asm.fa"
        infile.write_text(">contig1\nATCG\n")
        workdir = tmp_path / "wdir"
        workdir.mkdir()
        (tmp_path / "R1.fq").write_text("@r\nA\n+\nI\n")
        (tmp_path / "R2.fq").write_text("@r\nA\n+\nI\n")
        outfile = str(tmp_path / "polished.fasta")

        args = _make_args(
            tmp_path,
            method="pilon",
            infile=str(infile),
            workdir=str(workdir),
            iterations=5,
            outfile=outfile,
        )

        mock_result = MagicMock()
        mock_result.returncode = 0
        subprocess_calls = []

        def _fake_run(cmd, **kw):
            subprocess_calls.append(cmd)
            # Pilon writes polished1.fasta and an empty .changes file
            (workdir / "polished1.fasta").write_text(">contig1\nATCG\n")
            (workdir / "polished1.changes").write_text("")
            return mock_result

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", side_effect=_fake_run):
            with patch("AAFTF.polish.make_bwa_bam", return_value="asm.bwa.bam"):
                run(None, args)

        # Only one pilon invocation — converged at iteration 1
        pilon_calls = [c for c in subprocess_calls if c and c[0] == "pilon"]
        assert len(pilon_calls) == 1
        assert (tmp_path / "polished.fasta").exists()

    def test_runs_multiple_iterations_when_changes_found(self, tmp_path):
        """Pilon continues to the next iteration when changes are non-zero."""
        infile = tmp_path / "asm.fa"
        infile.write_text(">contig1\nATCG\n")
        workdir = tmp_path / "wdir"
        workdir.mkdir()
        (tmp_path / "R1.fq").write_text("@r\nA\n+\nI\n")
        (tmp_path / "R2.fq").write_text("@r\nA\n+\nI\n")
        outfile = str(tmp_path / "polished.fasta")

        args = _make_args(
            tmp_path,
            method="pilon",
            infile=str(infile),
            workdir=str(workdir),
            iterations=2,
            outfile=outfile,
        )

        mock_result = MagicMock()
        mock_result.returncode = 0
        call_count = [0]

        def _fake_run(cmd, **kw):
            call_count[0] += 1
            i = call_count[0]
            (workdir / f"polished{i}.fasta").write_text(">contig1\nATCG\n")
            # Non-empty changes on iteration 1, empty on iteration 2
            n_changes = 1 if i == 1 else 0
            (workdir / f"polished{i}.changes").write_text("change\n" * n_changes)
            return mock_result

        from AAFTF.polish import run

        with patch("AAFTF.polish.subprocess.run", side_effect=_fake_run):
            with patch("AAFTF.polish.make_bwa_bam", return_value="asm.bwa.bam"):
                run(None, args)

        assert call_count[0] == 2
        assert (tmp_path / "polished.fasta").exists()


# ---------------------------------------------------------------------------
# Integration test — real pilon run on Rhizopus test data
# ---------------------------------------------------------------------------

# Paths relative to the repository root
_TESTS_DIR = Path(__file__).parent
_INPUT_FASTA = _TESTS_DIR / "Rhizopus_microsporus_NRRL_5546.fcs_screen.fasta"
_R1 = _TESTS_DIR / "Rhizopus_microsporus_NRRL_5546_R1.fq.gz"
_R2 = _TESTS_DIR / "Rhizopus_microsporus_NRRL_5546_R2.fq.gz"

_have_pilon = shutil.which("pilon") is not None
_have_bwa = shutil.which("bwa") is not None
_test_data_present = _INPUT_FASTA.exists() and _R1.exists() and _R2.exists()

_skip_reason = []
if not _have_pilon:
    _skip_reason.append("pilon not in PATH")
if not _have_bwa:
    _skip_reason.append("bwa not in PATH")
if not _test_data_present:
    _skip_reason.append("test FASTA or reads missing")

_integration_skip = pytest.mark.skipif(
    bool(_skip_reason),
    reason=", ".join(_skip_reason) if _skip_reason else "",
)


@pytest.mark.integration
class TestPolishPilonIntegration:
    """Real pilon run using Rhizopus test data.

    Skipped automatically when pilon, bwa, or the test data files are absent.
    Run with:  pytest tests/test_polish.py -m integration -v
    """

    @_integration_skip
    def test_pilon_produces_output_fasta(self, tmp_path):
        """AAFTF polish --method pilon creates a non-empty polished FASTA."""
        outfile = str(tmp_path / "Rhizopus_microsporus_NRRL_5546.polish.fasta")
        args = Namespace(
            method="pilon",
            infile=str(_INPUT_FASTA),
            outfile=outfile,
            left=str(_R1),
            right=str(_R2),
            longreads=None,
            workdir=str(tmp_path / "workdir"),
            cpus=4,
            memory=16,
            iterations=1,
            diploid=False,
            ploidy=1,
            debug=False,
            pipe=True,
            polca="polca.sh",
        )

        from AAFTF.polish import run

        run(None, args)

        out = Path(outfile)
        assert out.exists(), "polished FASTA was not created"
        assert out.stat().st_size > 0, "polished FASTA is empty"

    @_integration_skip
    def test_pilon_output_is_valid_fasta(self, tmp_path):
        """Polished output contains at least one FASTA record."""
        outfile = str(tmp_path / "Rhizopus_microsporus_NRRL_5546.polish.fasta")
        args = Namespace(
            method="pilon",
            infile=str(_INPUT_FASTA),
            outfile=outfile,
            left=str(_R1),
            right=str(_R2),
            longreads=None,
            workdir=str(tmp_path / "workdir"),
            cpus=4,
            memory=16,
            iterations=1,
            diploid=False,
            ploidy=1,
            debug=False,
            pipe=True,
            polca="polca.sh",
        )

        from AAFTF.polish import run

        run(None, args)

        with open(outfile) as fh:
            first_line = fh.readline()
        assert first_line.startswith(">"), "output is not a FASTA file"

    @_integration_skip
    def test_pilon_output_size_reasonable(self, tmp_path):
        """Polished assembly should be within 10 % of input size."""
        outfile = str(tmp_path / "Rhizopus_microsporus_NRRL_5546.polish.fasta")
        args = Namespace(
            method="pilon",
            infile=str(_INPUT_FASTA),
            outfile=outfile,
            left=str(_R1),
            right=str(_R2),
            longreads=None,
            workdir=str(tmp_path / "workdir"),
            cpus=4,
            memory=16,
            iterations=1,
            diploid=False,
            ploidy=1,
            debug=False,
            pipe=True,
            polca="polca.sh",
        )

        from AAFTF.polish import run

        run(None, args)

        input_size = _INPUT_FASTA.stat().st_size
        output_size = Path(outfile).stat().st_size
        ratio = output_size / input_size
        assert 0.9 <= ratio <= 1.1, f"Output size {output_size:,} is unexpectedly far from input {input_size:,} (ratio {ratio:.2f})"
