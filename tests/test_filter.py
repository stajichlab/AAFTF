"""Unit tests for aaftf/filter.py.

Covers:
  - CLI parser defaults and flag presence for the 'filter' subcommand
  - run() guard: no read1 reads → ValueError
  - bbduk command construction for paired-end and single-end reads
  - bwa/bowtie2/minimap2 command construction
  - contamdb FASTA is created from source files
  - basename auto-derivation from read1-reads filename

External network access and actual aligners are mocked throughout.
"""

import gzip
import os
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.main import main
from aaftf.resources import SEQ_DBS
from aaftf.utility import db_write_dir

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF filter ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_filter(argv):
    captured = {}

    def _capture(**kwargs):
        captured["args"] = Namespace(**kwargs)

    with patch.object(sys, "argv", argv):
        with patch("aaftf.filter.run", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for filter.run()
# ---------------------------------------------------------------------------


_UNSET = object()


def _make_filter_args(tmp_path, read1=_UNSET, read2=None, aligner="bbduk", **overrides):
    if read1 is _UNSET:
        read1 = str(tmp_path / "sample_R1.fastq.gz")
    defaults = dict(
        workdir=str(tmp_path / "workdir"),
        cpus=1,
        memory=None,
        read1=read1,
        read2=read2,
        basename=None,
        aligner=aligner,
        screen_accessions=None,
        screen_urls=None,
        screen_local=None,
        debug=False,
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


def _mock_download(url, dest, force=False):
    """Side effect for download_file: create a stub file at *dest* and return it."""
    p = Path(dest)
    p.parent.mkdir(parents=True, exist_ok=True)
    if dest.endswith(".gz"):
        _write_stub_gz(p)
    else:
        _write_stub_plain(p)
    return dest


@pytest.fixture(autouse=True)
def _stub_databases():
    """Put a stub contamdb.fa in the (test-isolated) database folder, as `AAFTF database contamdb` would."""
    Path(db_write_dir(), "contamdb.fa").write_text(">univec_core\nACGT\n>phix\nGGCC\n>lambda\nTTGG\n")


@pytest.fixture(autouse=True)
def _stub_dedupe(monkeypatch):
    """dedupe.sh is not run in unit tests: the combined FASTA is copied through unchanged."""

    def _copy(src, dest, debug=False):
        Path(dest).write_text(Path(src).read_text())
        return Path(dest).read_text().count(">")

    monkeypatch.setattr("aaftf.filter.dedupe_fasta", _copy)


def _fake_bbtools(root: Path) -> Path:
    """Lay out a BBTools install (bin/bbduk.sh -> opt/bbmap/bbduk.sh, opt/bbmap/resources); return bin/."""
    opt = root / "opt" / "bbmap"
    (opt / "resources").mkdir(parents=True)
    (opt / "bbduk.sh").write_text("#!/bin/sh\n")
    (opt / "bbduk.sh").chmod(0o755)
    _write_stub_gz(opt / "resources" / "sequencing_artifacts.fa.gz", b">contam_0\nACGT\n")
    _write_stub_gz(opt / "resources" / "lambda.fa.gz", b">lambda\nTTGG\n")
    bindir = root / "bin"
    bindir.mkdir()
    (bindir / "bbduk.sh").symlink_to(opt / "bbduk.sh")
    return bindir


@pytest.fixture(autouse=True)
def _stub_bbtools(tmp_path_factory, monkeypatch):
    """Put a fake BBTools install (with its bundled artifact/lambda FASTAs) first on PATH."""
    bindir = _fake_bbtools(tmp_path_factory.mktemp("bbtools"))
    monkeypatch.setenv("PATH", f"{bindir}:{__import__('os').environ['PATH']}")


# ---------------------------------------------------------------------------
# Parser defaults and required flags
# ---------------------------------------------------------------------------


class TestFilterParser:
    def test_debug_false_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "-v"])
        assert args.debug is True

    def test_default_aligner_is_bbduk(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.aligner == "bbduk"

    def test_aligner_bowtie2(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "--aligner", "bowtie2"])
        assert args.aligner == "bowtie2"

    def test_aligner_bwa(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "--aligner", "bwa"])
        assert args.aligner == "bwa"

    def test_aligner_minimap2(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "--aligner", "minimap2"])
        assert args.aligner == "minimap2"

    def test_default_cpus(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "-c", "4"])
        assert args.cpus == 4

    def test_screen_accessions_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.screen_accessions is None

    def test_screen_urls_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.screen_urls is None

    def test_screen_local_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.screen_local is None

    def test_parses_read1_reads(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.read1 == "R1.fq"

    def test_parses_read2_reads(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "-2", "R2.fq"])
        assert args.read2 == "R2.fq"

    def test_read2_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.read2 is None

    def test_missing_read1_exits_nonzero(self):
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
# run() guard: no read1 reads
# ---------------------------------------------------------------------------


class TestFilterRunGuards:
    def test_no_read1_raises(self, tmp_path):
        args = _make_filter_args(tmp_path, read1=None)
        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=0):
                with pytest.raises(ValueError):
                    run(**vars(args))


# ---------------------------------------------------------------------------
# contamdb creation
# ---------------------------------------------------------------------------


class TestFilterContamdbCreation:
    def test_contamdb_created_in_workdir(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        args = _make_filter_args(tmp_path, read1=read1, aligner="bbduk")
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)

        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))

        contamdb = workdir / "contamdb.fa"
        assert contamdb.exists()

    def test_screen_local_added_to_contamdb(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        local_fa = tmp_path / "extra.fa"
        local_fa.write_text(">extra\nATCGATCG\n")
        args = _make_filter_args(tmp_path, read1=read1, aligner="bbduk", screen_local=[str(local_fa)])
        Path(args.workdir).mkdir(parents=True, exist_ok=True)

        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))

        contamdb = Path(args.workdir) / "contamdb.fa"
        content = contamdb.read_text()
        assert ">extra" in content
        assert not (Path(args.workdir) / "contamdb.raw.fa").exists()  # the combined intermediate is removed

    def test_contamdb_contains_downloaded_stubs(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        args = _make_filter_args(tmp_path, read1=read1, aligner="bbduk")
        Path(args.workdir).mkdir(parents=True, exist_ok=True)

        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))

        contamdb = Path(args.workdir) / "contamdb.fa"
        assert contamdb.stat().st_size > 0

    def test_contamdb_starts_from_the_database_fasta(self, tmp_path):
        """Every --aligner screens the shared contamdb database (UniVec_Core, PhiX, lambda, adapters)."""
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        args = _make_filter_args(tmp_path, read1=read1, aligner="bwa")
        Path(args.workdir).mkdir(parents=True, exist_ok=True)

        from aaftf.filter import run

        with patch("aaftf.filter.count_fastq", return_value=100), patch("aaftf.filter.run_cmd"), patch("aaftf.filter.align_to_sorted_bam"):
            with pytest.raises(RuntimeError):  # the mocked aligner writes no BAM
                run(**vars(args))

        content = (Path(args.workdir) / "contamdb.fa").read_text()
        assert content.startswith(">univec_core\nACGT\n>phix\nGGCC\n>lambda\n")


class TestBbtoolsResources:
    def test_found_next_to_resolved_bbduk(self):
        from aaftf.utility import bbtools_resources

        assert [Path(p).name for p in bbtools_resources(["sequencing_artifacts.fa.gz", "lambda.fa.gz"])] == ["sequencing_artifacts.fa.gz", "lambda.fa.gz"]
        assert all(Path(p).parent.name == "resources" for p in bbtools_resources(["lambda.fa.gz"]))

    def test_missing_bbduk_raises(self, monkeypatch):
        from aaftf.utility import bbtools_resources

        monkeypatch.setattr("aaftf.utility.shutil.which", lambda name: None)
        with pytest.raises(FileNotFoundError, match="bbduk.sh not found"):
            bbtools_resources(["lambda.fa.gz"])

    def test_missing_resource_raises(self, tmp_path, monkeypatch):
        from aaftf.utility import bbtools_resources

        bindir = _fake_bbtools(tmp_path)
        (tmp_path / "opt" / "bbmap" / "resources" / "lambda.fa.gz").unlink()
        monkeypatch.setattr("aaftf.utility.shutil.which", lambda name: str(bindir / name))
        with pytest.raises(FileNotFoundError, match="lambda.fa.gz"):
            bbtools_resources(["lambda.fa.gz"])


# ---------------------------------------------------------------------------
# bbduk command construction
# ---------------------------------------------------------------------------


def _run_filter_bbduk(tmp_path, read1, read2=None, **extra):
    """Run filter.run() with aligner=bbduk; return captured subprocess commands."""
    args = _make_filter_args(tmp_path, read1=read1, read2=read2, aligner="bbduk", **extra)
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    cmds = []

    from aaftf.filter import run

    with patch("aaftf.filter.download_file", side_effect=_mock_download):
        with patch("aaftf.filter.count_fastq", return_value=100):
            with patch("aaftf.utility.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                run(**vars(args))
    return cmds, args


class TestFilterRunBbduk:
    # PE bbduk runs as shuffle.sh -> bbduk.sh -> reformat.sh (see a136d93,
    # the BBDuk PairStreamer paired-mode workaround), so a given argument may
    # land on any of the three commands rather than the first.
    def test_pe_command_includes_in(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any(f"in1={read1}" in " ".join(c) for c in cmds)

    def test_pe_command_includes_in2(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any(f"in2={read2}" in " ".join(c) for c in cmds)

    def test_pe_output_files_use_filtered_basename(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any("sample_filtered_1.fastq.gz" in " ".join(c) for c in cmds)

    def test_se_output_uses_filtered_1(self, tmp_path):
        """Single-end output has the same name as read 1 of a pair, which the pipelines expect."""
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2=None)
        cmd_str = " ".join(cmds[0])
        assert "out=sample_filtered_1.fastq.gz" in cmd_str

    def test_se_hint_is_a_valid_assemble_command(self, tmp_path, caplog):
        """The single-end hint passes the reads as -1 (assemble requires --method and -1), not only --merged."""
        import logging

        with caplog.at_level(logging.INFO, logger="aaftf"):
            _run_filter_bbduk(tmp_path, str(tmp_path / "sample_R1.fastq.gz"), read2=None)
        hint = caplog.text.split("Your next command might be:")[-1]
        assert "AAFTF assemble --method spades -1 sample_filtered_1.fastq.gz -c" in hint and "--merged" not in hint

    def test_command_starts_with_bbduk(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any(c[0] == "bbduk.sh" for c in cmds)

    def test_command_includes_contamdb_ref(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, args = _run_filter_bbduk(tmp_path, read1, read2)
        # contamdb.fa (which already holds UniVec_Core, PhiX, lambda and adapters) is the only ref= of bbduk.sh
        refs = [arg for c in cmds if c[0] == "bbduk.sh" for arg in c if arg.startswith("ref=")]
        assert refs == [f"ref={Path(args.workdir, 'contamdb.fa')}"]

    def test_basename_derived_from_underscore_split(self, tmp_path):
        read1 = str(tmp_path / "MySample_R1.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1)
        assert any("MySample_filtered_1.fastq.gz" in " ".join(c) for c in cmds)

    def test_explicit_basename_preserved(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, basename="custom")
        assert any("custom_filtered_1.fastq.gz" in " ".join(c) for c in cmds)


# ---------------------------------------------------------------------------
# Aligner-based filtering (bwa / bowtie2 / minimap2)
# ---------------------------------------------------------------------------


def _run_filter_aligner(tmp_path, aligner, read1, read2=None, **extra):
    """Run filter.run() with an aligner; return (run_cmd commands, aligner commands, args).

    ``align_to_sorted_bam`` is faked: it records the aligner command and touches
    the BAM path so the post-alignment ``samtools fastq`` block runs.
    """
    args = _make_filter_args(tmp_path, read1=read1, read2=read2, aligner=aligner, **extra)
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    cmds = []
    align_cmds = []

    def _fake_align(cmd, bam, *a, **kw):
        align_cmds.append(cmd)
        Path(bam).touch()

    from aaftf.filter import run

    with patch("aaftf.filter.download_file", side_effect=_mock_download):
        with patch("aaftf.filter.count_fastq", return_value=100):
            with patch("aaftf.filter.run_cmd", side_effect=lambda cmd, *a, **kw: cmds.append(cmd)):
                with patch("aaftf.filter.align_to_sorted_bam", side_effect=_fake_align):
                    run(**vars(args))
    return cmds, align_cmds, args


_ALIGNERS = ["bwa", "bowtie2", "minimap2"]
_INDEX_CMD = {"bwa": ["bwa", "index"], "bowtie2": ["bowtie2-build"]}


class TestFilterRunAligners:
    @pytest.mark.parametrize("aligner", _ALIGNERS)
    def test_pe_samtools_fastq_extracts_unmapped_pairs(self, tmp_path, aligner):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _, args = _run_filter_aligner(tmp_path, aligner, read1, read2)
        bam = str(Path(args.workdir, "sample_contam_db.bam"))
        expected = ["samtools", "fastq", "-f", "12", "-1", "sample_filtered_1.fastq.gz", "-2", "sample_filtered_2.fastq.gz", bam]
        assert expected in cmds

    @pytest.mark.parametrize("aligner", _ALIGNERS)
    def test_se_samtools_fastq_extracts_unmapped_reads(self, tmp_path, aligner):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _, args = _run_filter_aligner(tmp_path, aligner, read1)
        bam = str(Path(args.workdir, "sample_contam_db.bam"))
        assert ["samtools", "fastq", "-f", "4", "-0", "sample_filtered_1.fastq.gz", bam] in cmds

    @pytest.mark.parametrize("aligner", ["bwa", "bowtie2"])
    def test_index_built(self, tmp_path, aligner):
        cmds, _, args = _run_filter_aligner(tmp_path, aligner, str(tmp_path / "sample_R1.fastq.gz"))
        contamdb = str(Path(args.workdir, "contamdb.fa"))
        idx = _INDEX_CMD[aligner]
        assert any(c[: len(idx)] == idx and contamdb in c for c in cmds)

    def test_bwa_command(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        _, align_cmds, _ = _run_filter_aligner(tmp_path, "bwa", read1, read2)
        assert align_cmds == [["bwa", "mem", "-t", "1", "contamdb.fa", read1, read2]]

    def test_bwa_se_command(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        _, align_cmds, _ = _run_filter_aligner(tmp_path, "bwa", read1)
        assert align_cmds == [["bwa", "mem", "-t", "1", "contamdb.fa", read1]]

    def test_bowtie2_pe_command(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        _, align_cmds, _ = _run_filter_aligner(tmp_path, "bowtie2", read1, read2)
        cmd = align_cmds[0]
        assert cmd[:3] == ["bowtie2", "-x", "contamdb.fa"]
        assert cmd[-4:] == ["-1", read1, "-2", read2]
        assert "-U" not in cmd

    def test_bowtie2_se_command(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        _, align_cmds, _ = _run_filter_aligner(tmp_path, "bowtie2", read1)
        cmd = align_cmds[0]
        assert cmd[-2:] == ["-U", read1]
        assert "-1" not in cmd

    def test_minimap2_pe_command(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        _, align_cmds, _ = _run_filter_aligner(tmp_path, "minimap2", read1, read2)
        assert align_cmds == [["minimap2", "-ax", "sr", "-t", "1", "contamdb.fa", read1, read2]]

    def test_minimap2_se_command(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        _, align_cmds, _ = _run_filter_aligner(tmp_path, "minimap2", read1)
        assert align_cmds == [["minimap2", "-ax", "sr", "-t", "1", "contamdb.fa", read1]]

    def test_rerun_reuses_contamdb_and_bam(self, tmp_path):
        """A rerun in the same work directory with unchanged sources rebuilds and realigns nothing."""
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        _run_filter_aligner(tmp_path, "bwa", read1)
        contamdb = tmp_path / "workdir" / "contamdb.fa"
        built = contamdb.stat().st_mtime_ns
        cmds, align_cmds, _ = _run_filter_aligner(tmp_path, "bwa", read1)
        assert align_cmds == []
        assert not any(c[:2] == ["bwa", "index"] for c in cmds)
        assert contamdb.stat().st_mtime_ns == built
        assert any(c[:2] == ["samtools", "fastq"] for c in cmds)

    def test_stale_bam_is_realigned(self, tmp_path):
        """A BAM older than contamdb.fa (e.g. left from an earlier database) is not reused."""
        workdir = tmp_path / "workdir"
        workdir.mkdir()
        bam = workdir / "sample_contam_db.bam"
        bam.touch()
        os.utime(bam, (0, 0))
        _, align_cmds, _ = _run_filter_aligner(tmp_path, "bwa", str(tmp_path / "sample_R1.fastq.gz"))
        assert len(align_cmds) == 1

    def test_changed_sources_rebuild_contamdb_and_realign(self, tmp_path):
        """Adding a screening source rebuilds contamdb.fa even though no source file is newer."""
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        extra = tmp_path / "extra.fa"
        extra.write_text(">extra\nACGT\n")
        os.utime(extra, (0, 0))
        _run_filter_aligner(tmp_path, "bwa", read1)
        _, align_cmds, args = _run_filter_aligner(tmp_path, "bwa", read1, screen_local=[str(extra)])
        assert ">extra" in Path(args.workdir, "contamdb.fa").read_text()
        assert len(align_cmds) == 1


class TestFilterDatabases:
    def _run(self, tmp_path, **extra):
        args = _make_filter_args(tmp_path, read1=str(tmp_path / "sample_R1.fastq.gz"), **extra)
        from aaftf.filter import run

        calls = []
        with patch("aaftf.filter.download_file", side_effect=lambda url, dest, force=False: (calls.append(url), _mock_download(url, dest))[1]):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))
        return calls

    def test_missing_contamdb_raises_with_build_command(self, tmp_path):
        Path(db_write_dir(), "contamdb.fa").unlink()
        with pytest.raises(FileNotFoundError) as exc:
            self._run(tmp_path)
        assert "AAFTF database contamdb" in str(exc.value)

    def test_databases_are_not_downloaded_by_filter(self, tmp_path):
        assert self._run(tmp_path) == []

    def test_screen_accession_fetched_from_eutils(self, tmp_path):
        calls = self._run(tmp_path, screen_accessions=["NC_001422"])
        assert calls == [SEQ_DBS["nucleotide"] % "NC_001422"]
        assert calls[0].startswith("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?")


def test_unknown_aligner_raises_before_any_work(tmp_path, monkeypatch):
    """An aligner outside the CLI choices is rejected up front (no work directory created)."""
    from aaftf.filter import run

    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="Unknown aligner"):
        run(read1=str(tmp_path / "R1.fq"), aligner="novoalign")
    assert not list(tmp_path.glob("aaftf-filter_*"))
