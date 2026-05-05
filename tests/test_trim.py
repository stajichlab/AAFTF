"""Unit tests for AAFTF/trim.py.

Covers:
  - CLI parser defaults and flag presence for the 'trim' subcommand
  - bbduk command construction for paired-end and single-end reads
  - fastp command construction including merge, dedup, and cut flags
  - trimmomatic guard: exits when no jar is found
  - basename auto-derivation from the left-reads filename
"""

import sys
from argparse import Namespace
from unittest.mock import patch

import pytest

from AAFTF.AAFTF_main import main

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF trim ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_trim(argv):
    captured = {}

    def _capture(parser, args):
        captured["args"] = args

    with patch.object(sys, "argv", argv):
        with patch("AAFTF.AAFTF_main.run_subtool", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for trim.run()
# ---------------------------------------------------------------------------


_UNSET = object()


def _make_trim_args(tmp_path, method="bbduk", left=None, right=_UNSET, **overrides):
    if left is None:
        left = str(tmp_path / "sample_R1.fastq.gz")
    if right is _UNSET:
        right = str(tmp_path / "sample_R2.fastq.gz")
    defaults = dict(
        method=method,
        left=left,
        right=right,
        basename=None,
        cpus=1,
        memory=None,
        minlen=75,
        avgqual=10,
        merge=False,
        dedup=False,
        cutfront=False,
        cuttail=False,
        cutright=False,
        debug=False,
        pipe=True,
        trimmomatic=None,
        trimmomatic_adaptors="TruSeq3-PE.fa",
        trimmomatic_clip="ILLUMINACLIP:%s:2:30:10",
        trimmomatic_leadingwindow=3,
        trimmomatic_trailingwindow=3,
        trimmomatic_slidingwindow="4:15",
        trimmomatic_quality="phred33",
    )
    defaults.update(overrides)
    return Namespace(**defaults)


# ---------------------------------------------------------------------------
# Parser defaults and required flags
# ---------------------------------------------------------------------------


class TestTrimParser:
    def test_debug_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "-v"])
        assert args.debug is True

    def test_pipe_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.pipe is False

    def test_pipe_flag_sets_true(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "--pipe"])
        assert args.pipe is True

    def test_default_method_is_bbduk(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.method == "bbduk"

    def test_method_trimmomatic(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "--method", "trimmomatic"])
        assert args.method == "trimmomatic"

    def test_method_fastp(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "--method", "fastp"])
        assert args.method == "fastp"

    def test_default_minlen(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.minlen == 75

    def test_custom_minlen(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "-ml", "50"])
        assert args.minlen == 50

    def test_default_avgqual(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.avgqual == 10

    def test_custom_avgqual(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "-aq", "20"])
        assert args.avgqual == 20

    def test_default_cpus(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "-c", "8"])
        assert args.cpus == 8

    def test_memory_none_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.memory is None

    def test_merge_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.merge is False

    def test_merge_flag_sets_true(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "--merge"])
        assert args.merge is True

    def test_dedup_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.dedup is False

    def test_dedup_flag_sets_true(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "--dedup"])
        assert args.dedup is True

    def test_cutfront_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.cutfront is False

    def test_cuttail_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.cuttail is False

    def test_cutright_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.cutright is False

    def test_parses_left_reads(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.left == "R1.fq"

    def test_parses_right_reads(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq", "-r", "R2.fq"])
        assert args.right == "R2.fq"

    def test_right_none_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.right is None

    def test_trim_reads_alias(self):
        args = _parse_trim(["AAFTF", "trim_reads", "-l", "R1.fq"])
        assert args.left == "R1.fq"

    def test_read_trim_alias(self):
        args = _parse_trim(["AAFTF", "read_trim", "-l", "R1.fq"])
        assert args.left == "R1.fq"

    def test_missing_left_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "trim"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_trim_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "trim", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0

    def test_default_trimmomatic_adaptors(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.trimmomatic_adaptors == "TruSeq3-PE.fa"

    def test_default_trimmomatic_quality(self):
        args = _parse_trim(["AAFTF", "trim", "-l", "R1.fq"])
        assert args.trimmomatic_quality == "phred33"


# ---------------------------------------------------------------------------
# bbduk command construction
# ---------------------------------------------------------------------------


def _run_bbduk(tmp_path, left, right=None, **extra):
    """Invoke trim.run() with method=bbduk, return captured subprocess commands."""
    args = _make_trim_args(tmp_path, method="bbduk", left=left, right=right, **extra)
    cmds = []
    with patch("AAFTF.trim.countfastq", return_value=100):
        with patch("AAFTF.trim.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)) as _:
            from AAFTF.trim import run

            run(None, args)
    return cmds, args


class TestTrimRunBbduk:
    def test_pe_command_includes_in1(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, left, right)
        assert any(f"in1={left}" in " ".join(c) for c in cmds)

    def test_pe_command_includes_in2(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, left, right)
        assert any(f"in2={right}" in " ".join(c) for c in cmds)

    def test_pe_command_includes_out1_out2(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, args = _run_bbduk(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        assert f"out1={args.basename}_1P.fastq.gz" in cmd_str
        assert f"out2={args.basename}_2P.fastq.gz" in cmd_str

    def test_se_command_includes_in(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, left, right=None)
        assert any(f"in={left}" in " ".join(c) for c in cmds)

    def test_se_command_includes_out(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        cmds, args = _run_bbduk(tmp_path, left, right=None)
        assert any(f"out={args.basename}_1U.fastq.gz" in " ".join(c) for c in cmds)

    def test_command_includes_minlen(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, left, right, minlen=50)
        assert any("minlen=50" in " ".join(c) for c in cmds)

    def test_command_includes_avgqual(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        right = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, left, right, avgqual=20)
        assert any("maq=20" in " ".join(c) for c in cmds)

    def test_basename_derived_from_underscore_split(self, tmp_path):
        left = str(tmp_path / "MySample_R1.fastq.gz")
        _, args = _run_bbduk(tmp_path, left)
        assert args.basename == "MySample"

    def test_basename_derived_from_dot_split(self, tmp_path):
        left = str(tmp_path / "MySample.R1.fastq.gz")
        _, args = _run_bbduk(tmp_path, left)
        assert args.basename == "MySample"

    def test_explicit_basename_not_overridden(self, tmp_path):
        left = str(tmp_path / "sample_R1.fastq.gz")
        _, args = _run_bbduk(tmp_path, left, basename="mybase")
        assert args.basename == "mybase"


# ---------------------------------------------------------------------------
# fastp command construction
# ---------------------------------------------------------------------------


def _run_fastp(tmp_path, left, right=None, **extra):
    args = _make_trim_args(tmp_path, method="fastp", left=left, right=right, **extra)
    cmds = []
    with patch("AAFTF.trim.countfastq", return_value=100):
        with patch("AAFTF.trim.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
            from AAFTF.trim import run

            run(None, args)
    return cmds, args


class TestTrimRunFastp:
    def test_pe_command_includes_in1_in2(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        assert f"--in1={left}" in cmd_str
        assert f"--in2={right}" in cmd_str

    def test_pe_command_includes_out1_out2(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, args = _run_fastp(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        assert f"--out1={args.basename}_1P.fastq.gz" in cmd_str
        assert f"--out2={args.basename}_2P.fastq.gz" in cmd_str

    def test_merge_adds_merge_flag_and_output(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, args = _run_fastp(tmp_path, left, right, merge=True)
        cmd_str = " ".join(cmds[0])
        assert "--merge" in cmd_str
        assert f"--merged_out={args.basename}_MG.fastq.gz" in cmd_str

    def test_dedup_adds_dedup_flag(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, left, right, dedup=True)
        assert "--dedup" in cmds[0]

    def test_cutfront_adds_flag(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, left, right, cutfront=True)
        assert "--cut_front" in cmds[0]

    def test_cuttail_adds_flag(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, left, right, cuttail=True)
        assert "--cut_tail" in cmds[0]

    def test_cutright_adds_flag(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, left, right, cutright=True)
        assert "--cut_right" in cmds[0]

    def test_se_uses_in_out(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        cmds, args = _run_fastp(tmp_path, left, right=None)
        cmd_str = " ".join(cmds[0])
        assert f"--in={left}" in cmd_str

    def test_merge_not_added_when_false(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, left, right, merge=False)
        assert "--merge" not in cmds[0]

    def test_command_includes_minlen(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, left, right, minlen=50)
        assert "50" in cmds[0]

    def test_html_json_reports_included(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        cmds, args = _run_fastp(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        assert ".fastp.html" in cmd_str
        assert ".fastp.json" in cmd_str


# ---------------------------------------------------------------------------
# trimmomatic guard
# ---------------------------------------------------------------------------


class TestTrimRunTrimmomatic:
    def test_no_jar_exits(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        args = _make_trim_args(tmp_path, method="trimmomatic", left=left, right=right, trimmomatic=None)

        from AAFTF.trim import run

        with patch("AAFTF.trim.find_trimmomatic", return_value=False):
            with patch("AAFTF.trim.countfastq", return_value=100):
                with pytest.raises(SystemExit):
                    run(None, args)

    def test_jar_found_builds_pe_command(self, tmp_path):
        left = str(tmp_path / "s_R1.fastq.gz")
        right = str(tmp_path / "s_R2.fastq.gz")
        fake_jar = str(tmp_path / "trimmomatic.jar")
        fake_adaptor = str(tmp_path / "TruSeq3-PE.fa")
        (tmp_path / "TruSeq3-PE.fa").write_text(">adapt\nATCG\n")
        args = _make_trim_args(
            tmp_path,
            method="trimmomatic",
            left=left,
            right=right,
            trimmomatic_adaptors=fake_adaptor,
        )
        cmds = []
        with patch("AAFTF.trim.find_trimmomatic", return_value=fake_jar):
            with patch("AAFTF.trim.countfastq", return_value=100):
                with patch("AAFTF.trim.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                    with patch("AAFTF.trim.Fzip_inplace"):
                        with patch("AAFTF.trim.SafeRemove"):
                            from AAFTF.trim import run

                            run(None, args)
        assert len(cmds) > 0
        assert "PE" in cmds[0]
        assert fake_jar in cmds[0]
