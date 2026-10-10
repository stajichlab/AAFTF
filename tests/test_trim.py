"""Unit tests for aaftf/trim.py.

Covers:
  - CLI parser defaults and flag presence for the 'trim' subcommand
  - bbduk command construction for paired-end and single-end reads
  - fastp command construction including merge, dedup, and cut flags
  - trimmomatic: FileNotFoundError when no jar is found, ValueError without reads,
    PE/SE command construction and _find_trimmomatic jar lookup
  - basename auto-derivation from the read1-reads filename
"""

import sys
from argparse import Namespace
from unittest.mock import patch

import pytest

from aaftf.main import main

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF trim ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_trim(argv):
    captured = {}

    def _capture(**kwargs):
        captured["args"] = Namespace(**kwargs)

    with patch.object(sys, "argv", argv):
        with patch("aaftf.trim.run", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for trim.run()
# ---------------------------------------------------------------------------


_UNSET = object()


def _make_trim_args(tmp_path, method="bbduk", read1=None, read2=_UNSET, **overrides):
    if read1 is None:
        read1 = str(tmp_path / "sample_R1.fastq.gz")
    if read2 is _UNSET:
        read2 = str(tmp_path / "sample_R2.fastq.gz")
    defaults = dict(
        method=method,
        read1=read1,
        read2=read2,
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
        trimmomatic_adaptors="TruSeq3-PE.fa",
        trimmomatic_clip="2:30:10",
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
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "-v"])
        assert args.debug is True

    def test_default_method_is_bbduk(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.method == "bbduk"

    def test_method_trimmomatic(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "--method", "trimmomatic"])
        assert args.method == "trimmomatic"

    def test_method_fastp(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "--method", "fastp"])
        assert args.method == "fastp"

    def test_default_minlen(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.minlen == 75

    def test_custom_minlen(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "-ml", "50"])
        assert args.minlen == 50

    def test_default_avgqual(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.avgqual == 10

    def test_custom_avgqual(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "-aq", "20"])
        assert args.avgqual == 20

    def test_default_cpus(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "-c", "8"])
        assert args.cpus == 8

    def test_memory_default_is_8gb_int(self):
        # an int, so bbduk gets -Xmx8g (Java rejects -Xmx8.0g and -XmxNoneg)
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.memory == 8 and isinstance(args.memory, int)

    def test_merge_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.merge is False

    def test_merge_flag_sets_true(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "--merge"])
        assert args.merge is True

    def test_dedup_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.dedup is False

    def test_dedup_flag_sets_true(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "--dedup"])
        assert args.dedup is True

    def test_cutfront_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.cutfront is False

    def test_cuttail_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.cuttail is False

    def test_cutright_false_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.cutright is False

    def test_parses_read1_reads(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.read1 == "R1.fq"

    def test_parses_read2_reads(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq", "-2", "R2.fq"])
        assert args.read2 == "R2.fq"

    def test_read2_none_by_default(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.read2 is None

    def test_missing_read1_exits_nonzero(self):
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
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.trimmomatic_adaptors == "TruSeq3-PE.fa"

    def test_default_trimmomatic_quality(self):
        args = _parse_trim(["AAFTF", "trim", "-1", "R1.fq"])
        assert args.trimmomatic_quality == "phred33"


# ---------------------------------------------------------------------------
# bbduk command construction
# ---------------------------------------------------------------------------


def _run_bbduk(tmp_path, read1, read2=None, **extra):
    """Invoke trim.run() with method=bbduk, return captured subprocess commands."""
    args = _make_trim_args(tmp_path, method="bbduk", read1=read1, read2=read2, **extra)
    cmds = []
    with patch("aaftf.trim.count_fastq", return_value=100):
        with patch("aaftf.utility.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)) as _:
            from aaftf.trim import run

            run(**vars(args))
    return cmds, args


class TestTrimRunBbduk:
    def test_pe_command_includes_in1(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, read2)
        assert any(f"in1={read1}" in " ".join(c) for c in cmds)

    def test_pe_command_includes_in2(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, read2)
        assert any(f"in2={read2}" in " ".join(c) for c in cmds)

    def test_pe_command_includes_out1_out2(self, tmp_path):
        # PE bbduk runs as shuffle.sh -> bbduk.sh -> reformat.sh (see a136d93);
        # the final out1/out2 filenames are written by the reformat.sh step,
        # not the first command.
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, read2)
        assert any("out1=sample_1P.fastq.gz" in " ".join(c) for c in cmds)
        assert any("out2=sample_2P.fastq.gz" in " ".join(c) for c in cmds)

    def test_se_command_includes_in(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, read2=None)
        assert any(f"in={read1}" in " ".join(c) for c in cmds)

    def test_se_command_includes_out(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, read2=None)
        assert any("out=sample_1U.fastq.gz" in " ".join(c) for c in cmds)

    def test_command_includes_minlen(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, read2, minlen=50)
        assert any("minlen=50" in " ".join(c) for c in cmds)

    def test_command_includes_avgqual(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, read2, avgqual=20)
        assert any("maq=20" in " ".join(c) for c in cmds)

    def test_basename_derived_from_underscore_split(self, tmp_path):
        read1 = str(tmp_path / "MySample_R1.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1)
        assert any("out=MySample_1U.fastq.gz" in " ".join(c) for c in cmds)

    def test_basename_derived_from_dot_split(self, tmp_path):
        read1 = str(tmp_path / "MySample.R1.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1)
        assert any("out=MySample_1U.fastq.gz" in " ".join(c) for c in cmds)

    def test_explicit_basename_not_overridden(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_bbduk(tmp_path, read1, basename="mybase")
        assert any("out=mybase_1U.fastq.gz" in " ".join(c) for c in cmds)


# ---------------------------------------------------------------------------
# fastp command construction
# ---------------------------------------------------------------------------


def _run_fastp(tmp_path, read1, read2=None, **extra):
    args = _make_trim_args(tmp_path, method="fastp", read1=read1, read2=read2, **extra)
    cmds = []
    with patch("aaftf.trim.count_fastq", return_value=100):
        with patch("aaftf.utility.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
            from aaftf.trim import run

            run(**vars(args))
    return cmds, args


class TestTrimRunFastp:
    def test_pe_command_includes_in1_in2(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2)
        cmd_str = " ".join(cmds[0])
        assert f"--in1={read1}" in cmd_str
        assert f"--in2={read2}" in cmd_str

    def test_pe_command_includes_out1_out2(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2)
        cmd_str = " ".join(cmds[0])
        assert "--out1=s_1P.fastq.gz" in cmd_str
        assert "--out2=s_2P.fastq.gz" in cmd_str

    def test_merge_adds_merge_flag_and_output(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2, merge=True)
        cmd_str = " ".join(cmds[0])
        assert "--merge" in cmd_str
        assert "--merged_out=s_MG.fastq.gz" in cmd_str

    def test_dedup_adds_dedup_flag(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2, dedup=True)
        assert "--dedup" in cmds[0]

    def test_cutfront_adds_flag(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2, cutfront=True)
        assert "--cut_front" in cmds[0]

    def test_cuttail_adds_flag(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2, cuttail=True)
        assert "--cut_tail" in cmds[0]

    def test_cutright_adds_flag(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2, cutright=True)
        assert "--cut_right" in cmds[0]

    def test_se_uses_in_out(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        cmds, args = _run_fastp(tmp_path, read1, read2=None)
        cmd_str = " ".join(cmds[0])
        assert f"--in={read1}" in cmd_str

    def test_merge_not_added_when_false(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2, merge=False)
        assert "--merge" not in cmds[0]

    def test_command_includes_minlen(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, _ = _run_fastp(tmp_path, read1, read2, minlen=50)
        assert "50" in cmds[0]

    def test_html_json_reports_included(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        cmds, args = _run_fastp(tmp_path, read1, read2)
        cmd_str = " ".join(cmds[0])
        assert ".fastp.html" in cmd_str
        assert ".fastp.json" in cmd_str


# ---------------------------------------------------------------------------
# trimmomatic guard
# ---------------------------------------------------------------------------


class TestTrimRunTrimmomatic:
    def test_no_jar_raises(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        args = _make_trim_args(tmp_path, method="trimmomatic", read1=read1, read2=read2)

        from aaftf.trim import run

        with patch("aaftf.trim._find_trimmomatic", return_value=None):
            with patch("aaftf.trim.count_fastq", return_value=100):
                with pytest.raises(FileNotFoundError):
                    run(**vars(args))

    def test_jar_found_builds_pe_command(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        read2 = str(tmp_path / "s_R2.fastq.gz")
        fake_jar = str(tmp_path / "trimmomatic.jar")
        fake_adaptor = str(tmp_path / "TruSeq3-PE.fa")
        (tmp_path / "TruSeq3-PE.fa").write_text(">adapt\nATCG\n")
        args = _make_trim_args(
            tmp_path,
            method="trimmomatic",
            read1=read1,
            read2=read2,
            trimmomatic_adaptors=fake_adaptor,
        )
        cmds = []
        with patch("aaftf.trim._find_trimmomatic", return_value=fake_jar):
            with patch("aaftf.trim.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                    with patch("aaftf.trim.safe_remove"):
                        from aaftf.trim import run

                        run(**vars(args))
        assert len(cmds) > 0
        assert "PE" in cmds[0]
        assert fake_jar in cmds[0]
        # trimmomatic writes gzipped output directly (by the .gz extension)
        assert any(a.endswith("_1P.fastq.gz") for a in cmds[0])
        assert any(a.endswith("_2P.fastq.gz") for a in cmds[0])

    def test_jar_found_builds_se_command(self, tmp_path):
        read1 = str(tmp_path / "s_R1.fastq.gz")
        fake_jar = str(tmp_path / "trimmomatic.jar")
        fake_adaptor = str(tmp_path / "TruSeq3-SE.fa")
        (tmp_path / "TruSeq3-SE.fa").write_text(">adapt\nATCG\n")
        args = _make_trim_args(tmp_path, method="trimmomatic", read1=read1, read2=None, basename="s", trimmomatic_adaptors=fake_adaptor)
        cmds = []
        with patch("aaftf.trim._find_trimmomatic", return_value=fake_jar):
            with patch("aaftf.trim.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                    from aaftf.trim import run

                    run(**vars(args))
        cmd = cmds[0]
        assert cmd[3] == "SE"
        assert "PE" not in cmd
        assert cmd[cmd.index(read1) + 1] == "s_1U.fastq.gz"

    def test_no_reads_raises(self, tmp_path):
        from aaftf.trim import run_trimmomatic

        adaptor = tmp_path / "TruSeq3-PE.fa"
        adaptor.write_text(">adapt\nATCG\n")
        with patch("aaftf.trim._find_trimmomatic", return_value=str(tmp_path / "trimmomatic.jar")) as find_jar:
            with patch("aaftf.trim.run_cmd") as run_cmd:
                with pytest.raises(ValueError, match="Must provide read1"):
                    run_trimmomatic(None, None, "s", 1, 75, str(adaptor), "2:30:10", 3, 3, "4:15", "phred33", False, False)
        find_jar.assert_not_called()  # checked before looking up the jar/adaptors
        run_cmd.assert_not_called()


class TestFindTrimmomatic:
    def _find(self, tmp_path, content):
        script = tmp_path / "bin" / "trimmomatic"
        script.parent.mkdir()
        script.write_text(content)
        from aaftf.trim import _find_trimmomatic

        with patch("aaftf.trim.shutil.which", return_value=str(script)):
            return _find_trimmomatic(), script

    def test_bash_wrapper_returns_jar_from_exec_line(self, tmp_path):
        found, _ = self._find(tmp_path, '#!/bin/bash\nset -e\nexec java -Xmx1g -jar /x/trimmomatic.jar "$@"\n')
        assert found == "/x/trimmomatic.jar"

    def test_python_wrapper_returns_jar_next_to_script(self, tmp_path):
        found, script = self._find(tmp_path, "#!/usr/bin/env python\nprint('hi')\n")
        assert found == str(script.resolve().parent / "trimmomatic.jar")

    def test_unknown_shebang_returns_none(self, tmp_path):
        found, _ = self._find(tmp_path, "#!/bin/sh\nexec java -jar /x/trimmomatic.jar\n")
        assert found is None

    def test_not_on_path_returns_none(self):
        from aaftf.trim import _find_trimmomatic

        with patch("aaftf.trim.shutil.which", return_value=None):
            assert _find_trimmomatic() is None


class TestTrimmomaticAdaptorLookup:
    """run_trimmomatic finds the adaptors next to the jar or under <prefix>/share/trimmomatic."""

    def _run(self, tmp_path, jar):
        args = _make_trim_args(tmp_path, method="trimmomatic", read1=str(tmp_path / "s_R1.fastq.gz"), read2=str(tmp_path / "s_R2.fastq.gz"), trimmomatic_adaptors="missing.fa")
        cmds = []
        with patch("aaftf.trim._find_trimmomatic", return_value=str(jar)):
            with patch("aaftf.trim.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                    with patch("aaftf.trim.safe_remove"):
                        from aaftf.trim import run

                        run(**vars(args))
        return next(a for a in cmds[0] if a.startswith("ILLUMINACLIP:"))

    def _adaptor(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(">adapt\nATCG\n")
        return path

    def test_adaptors_next_to_jar(self, tmp_path):
        adaptors = self._adaptor(tmp_path / "share" / "trimmomatic" / "adapters" / "TruSeq3-PE.fa")
        assert self._run(tmp_path, tmp_path / "share" / "trimmomatic" / "trimmomatic.jar").startswith(f"ILLUMINACLIP:{adaptors}:")

    def test_adaptors_under_prefix_share(self, tmp_path):
        adaptors = self._adaptor(tmp_path / "prefix" / "share" / "trimmomatic" / "adapters" / "TruSeq3-PE.fa")
        assert self._run(tmp_path, tmp_path / "prefix" / "lib" / "java" / "trimmomatic.jar").startswith(f"ILLUMINACLIP:{adaptors}:")

    def test_missing_adaptors_raise(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="--trimmomatic_adaptors"):
            self._run(tmp_path, tmp_path / "lib" / "trimmomatic.jar")
