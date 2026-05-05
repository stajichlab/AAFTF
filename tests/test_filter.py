"""Unit tests for AAFTF/filter.py.

Covers:
  - CLI parser defaults and flag presence for the 'filter' subcommand
  - run() guard: no left reads → sys.exit(1)
  - bbduk command construction for paired-end and single-end reads
  - bwa/bowtie2/minimap2 command construction
  - contamdb FASTA is created from source files
  - basename auto-derivation from left-reads filename

External network access and actual aligners are mocked throughout.
"""

import gzip
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from AAFTF.AAFTF_main import main

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF filter ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_filter(argv):
    captured = {}

    def _capture(parser, args):
        captured["args"] = args

    with patch.object(sys, "argv", argv):
        with patch("AAFTF.AAFTF_main.run_subtool", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for filter.run()
# ---------------------------------------------------------------------------


_UNSET = object()


def _make_filter_args(tmp_path, left=_UNSET, right=None, aligner="bbduk", **overrides):
    if left is _UNSET:
        left = str(tmp_path / "sample_R1.fastq.gz")
    defaults = dict(
        workdir=str(tmp_path / "workdir"),
        cpus=1,
        memory=None,
        left=left,
        right=right,
        basename=None,
        AAFTF_DB=None,
        aligner=aligner,
        screen_accessions=None,
        screen_urls=None,
        screen_local=None,
        debug=False,
        pipe=True,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


# ---------------------------------------------------------------------------
# Create stub contaminant files (valid gzip + plain) so contamdb can be built
# ---------------------------------------------------------------------------


def _write_stub_gz(path: Path, content: bytes = b">stub\nATCG\n"):
    with gzip.open(path, "wb") as f:
        f.write(content)


def _write_stub_plain(path: Path, content: bytes = b">stub\nATCG\n"):
    path.write_bytes(content)


def _mock_urlretrieve(url, dest):
    """Side effect for urlretrieve: create a stub file at *dest*."""
    p = Path(dest)
    p.parent.mkdir(parents=True, exist_ok=True)
    if dest.endswith(".gz"):
        _write_stub_gz(p)
    else:
        _write_stub_plain(p)


# ---------------------------------------------------------------------------
# Parser defaults and required flags
# ---------------------------------------------------------------------------


class TestFilterParser:
    def test_debug_false_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq", "-v"])
        assert args.debug is True

    def test_pipe_false_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.pipe is False

    def test_pipe_flag_sets_true(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq", "--pipe"])
        assert args.pipe is True

    def test_default_aligner_is_bbduk(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.aligner == "bbduk"

    def test_aligner_bowtie2(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq", "--aligner", "bowtie2"])
        assert args.aligner == "bowtie2"

    def test_aligner_bwa(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq", "--aligner", "bwa"])
        assert args.aligner == "bwa"

    def test_aligner_minimap2(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq", "--aligner", "minimap2"])
        assert args.aligner == "minimap2"

    def test_default_cpus(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq", "-c", "4"])
        assert args.cpus == 4

    def test_screen_accessions_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.screen_accessions is None

    def test_screen_urls_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.screen_urls is None

    def test_screen_local_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.screen_local is None

    def test_AAFTF_DB_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.AAFTF_DB is None

    def test_parses_left_reads(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.left == "R1.fq"

    def test_parses_right_reads(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq", "-r", "R2.fq"])
        assert args.right == "R2.fq"

    def test_right_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-l", "R1.fq"])
        assert args.right is None

    def test_filter_reads_alias(self):
        args = _parse_filter(["AAFTF", "filter_reads", "-l", "R1.fq"])
        assert args.left == "R1.fq"

    def test_read_filter_alias(self):
        args = _parse_filter(["AAFTF", "read_filter", "-l", "R1.fq"])
        assert args.left == "R1.fq"

    def test_missing_left_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "filter"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_filter_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "filter", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0


# ---------------------------------------------------------------------------
# run() guard: no left reads
# ---------------------------------------------------------------------------


class TestFilterRunGuards:
    def test_no_left_exits(self, tmp_path):
        args = _make_filter_args(tmp_path, left=None)
        from AAFTF.filter import run

        with patch("AAFTF.filter.urllib.request.urlretrieve", side_effect=_mock_urlretrieve):
            with patch("AAFTF.filter.countfastq", return_value=0):
                with pytest.raises(SystemExit):
                    run(None, args)


# ---------------------------------------------------------------------------
# contamdb creation
# ---------------------------------------------------------------------------


class TestFilterContamdbCreation:
    def test_contamdb_created_in_workdir(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        args = _make_filter_args(tmp_path, left=left, aligner="bbduk")
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)

        from AAFTF.filter import run

        with patch("AAFTF.filter.urllib.request.urlretrieve", side_effect=_mock_urlretrieve):
            with patch("AAFTF.filter.countfastq", return_value=100):
                with patch("AAFTF.filter.subprocess.run"):
                    run(None, args)

        contamdb = workdir / "contamdb.fa"
        assert contamdb.exists()

    def test_screen_local_added_to_contamdb(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        local_fa = tmp_path / "extra.fa"
        local_fa.write_text(">extra\nATCGATCG\n")
        args = _make_filter_args(tmp_path, left=left, aligner="bbduk", screen_local=[str(local_fa)])
        Path(args.workdir).mkdir(parents=True, exist_ok=True)

        from AAFTF.filter import run

        with patch("AAFTF.filter.urllib.request.urlretrieve", side_effect=_mock_urlretrieve):
            with patch("AAFTF.filter.countfastq", return_value=100):
                with patch("AAFTF.filter.subprocess.run"):
                    run(None, args)

        contamdb = Path(args.workdir) / "contamdb.fa"
        content = contamdb.read_text()
        assert ">extra" in content

    def test_contamdb_contains_downloaded_stubs(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        args = _make_filter_args(tmp_path, left=left, aligner="bbduk")
        Path(args.workdir).mkdir(parents=True, exist_ok=True)

        from AAFTF.filter import run

        with patch("AAFTF.filter.urllib.request.urlretrieve", side_effect=_mock_urlretrieve):
            with patch("AAFTF.filter.countfastq", return_value=100):
                with patch("AAFTF.filter.subprocess.run"):
                    run(None, args)

        contamdb = Path(args.workdir) / "contamdb.fa"
        assert contamdb.stat().st_size > 0


# ---------------------------------------------------------------------------
# bbduk command construction
# ---------------------------------------------------------------------------


def _run_filter_bbduk(tmp_path, left, right=None, **extra):
    """Run filter.run() with aligner=bbduk; return captured subprocess commands."""
    args = _make_filter_args(tmp_path, left=left, right=right, aligner="bbduk", **extra)
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    cmds = []

    from AAFTF.filter import run

    with patch("AAFTF.filter.urllib.request.urlretrieve", side_effect=_mock_urlretrieve):
        with patch("AAFTF.filter.countfastq", return_value=100):
            with patch("AAFTF.filter.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                run(None, args)
    return cmds, args


class TestFilterRunBbduk:
    def test_pe_command_includes_in(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, left, right)
        bbduk_cmd = cmds[0]
        cmd_str = " ".join(bbduk_cmd)
        assert f"in={left}" in cmd_str

    def test_pe_command_includes_in2(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        assert f"in2={right}" in cmd_str

    def test_pe_output_files_use_filtered_basename(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, args = _run_filter_bbduk(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        expected_out1 = f"{args.basename}_filtered_1.fastq.gz"
        assert expected_out1 in cmd_str

    def test_se_output_uses_U_suffix(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        cmds, args = _run_filter_bbduk(tmp_path, left, right=None)
        cmd_str = " ".join(cmds[0])
        assert f"{args.basename}_filtered_U.fastq.gz" in cmd_str

    def test_command_starts_with_bbduk(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, left, right)
        assert cmds[0][0] == "bbduk.sh"

    def test_command_includes_contamdb_ref(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, args = _run_filter_bbduk(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        # contamdb.fa is passed via ref= argument
        assert "contamdb.fa" in cmd_str

    def test_basename_derived_from_underscore_split(self, tmp_path):
        left = str(tmp_path / "MySample_R1.fastq.gz")
        _, args = _run_filter_bbduk(tmp_path, left)
        assert args.basename == "MySample"

    def test_explicit_basename_preserved(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        _, args = _run_filter_bbduk(tmp_path, left, basename="custom")
        assert args.basename == "custom"


# ---------------------------------------------------------------------------
# bwa command construction
# ---------------------------------------------------------------------------


def _run_filter_bwa(tmp_path, left, right=None, **extra):
    """Run filter.run() with aligner=bwa; return captured subprocess commands."""
    args = _make_filter_args(tmp_path, left=left, right=right, aligner="bwa", **extra)
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    cmds = []
    popen_cmds = []

    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.communicate.return_value = (b"", b"")

    def _fake_run(cmd, **kw):
        cmds.append(cmd)
        # When samtools sort is called, create the alignBAM so post-processing triggers
        if cmd and len(cmd) > 1 and "samtools" in str(cmd[0]):
            # find the output BAM argument (samtools sort -o <out>)
            for i, tok in enumerate(cmd):
                if tok == "-o" and i + 1 < len(cmd):
                    Path(cmd[i + 1]).touch()
                    break
            # fallback: touch any .bam in workdir
            for f in workdir.glob("*.bam"):
                pass  # already exists via above
            else:
                # create a stub alignBAM so isfile check passes
                (workdir / "_stub.bam").touch()

    from packaging.version import Version

    from AAFTF.filter import run

    with patch("AAFTF.filter.urllib.request.urlretrieve", side_effect=_mock_urlretrieve):
        with patch("AAFTF.filter.countfastq", return_value=100):
            with patch("AAFTF.filter.subprocess.run", side_effect=_fake_run):
                with patch("AAFTF.filter.subprocess.Popen", side_effect=lambda cmd, **kw: (popen_cmds.append(cmd), mock_proc)[1]):
                    with patch("AAFTF.utility.get_samtools_version", return_value=Version("1.23")):
                        with patch("AAFTF.filter.bam_read_count", return_value=(50, 50)):
                            with patch("AAFTF.filter.SafeRemove"):
                                run(None, args)
    return cmds, popen_cmds, args


class TestFilterRunBwa:
    def test_bwa_index_called(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, popen_cmds, _ = _run_filter_bwa(tmp_path, left, right)
        all_cmds = cmds + popen_cmds
        assert any(c and c[0] == "bwa" and "index" in c for c in all_cmds)

    def test_bwa_mem_called(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, popen_cmds, _ = _run_filter_bwa(tmp_path, left, right)
        all_cmds = cmds + popen_cmds
        assert any(c and c[0] == "bwa" and "mem" in c for c in all_cmds)

    def test_bwa_mem_includes_reads(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, popen_cmds, _ = _run_filter_bwa(tmp_path, left, right)
        mem_cmds = [c for c in popen_cmds if c and c[0] == "bwa" and "mem" in c]
        assert len(mem_cmds) > 0
        assert left in mem_cmds[0]
