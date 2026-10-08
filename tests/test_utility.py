"""Unit tests for aaftf/utility.py.

All tests are pure Python — no external bioinformatics tools required.
"""

import gzip
import io
import os
import shutil
import subprocess
from unittest.mock import patch

import pytest

from aaftf.utility import (
    COMPLEMENT,
    PafHit,
    align_to_sorted_bam,
    available_cpus,
    bam_read_count,
    basename_from_reads,
    calc_nx,
    check_file,
    cleanup_workdir,
    concat_files,
    count_fastq,
    db_dirs,
    db_file,
    db_write_dir,
    download_file,
    execute,
    fasta_stats,
    filter_fasta,
    find_db_file,
    find_db_index,
    home_db_cache,
    make_workdir,
    next_step_name,
    open_maybe_gz,
    paf_hits,
    require_databases,
    require_tools,
    run_cmd,
    safe_remove,
    samtools_sort_cmd,
    setup_logging,
    softwrap,
    warn_if_home_cache,
    write_fasta,
)
from tests.conftest import make_fastq_text

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# calc_nx
# ---------------------------------------------------------------------------


class TestCalcNx:
    # lengths [20, 80, 150] → total=250; 50% = 125 → 150 alone covers it
    def test_n50_l50(self):
        assert calc_nx([20, 80, 150]) == (150, 1)

    # 90% = 225 → 150 < 225, 150+80 = 230 ≥ 225
    def test_n90_l90(self):
        assert calc_nx([20, 80, 150], 0.9) == (80, 2)

    def test_single_contig(self):
        assert calc_nx([500]) == (500, 1)

    def test_equal_contigs(self):
        # [100]*4 → total=400, 50% = 200 → reached at the 2nd contig
        assert calc_nx([100, 100, 100, 100]) == (100, 2)

    def test_input_order_irrelevant(self):
        assert calc_nx([150, 20, 80]) == calc_nx([20, 80, 150])

    def test_input_not_modified(self):
        lengths = [20, 150, 80]
        calc_nx(lengths)
        assert lengths == [20, 150, 80]

    def test_empty(self):
        assert calc_nx([]) == (0, 0)


# ---------------------------------------------------------------------------
# softwrap
# ---------------------------------------------------------------------------


class TestSoftwrap:
    def test_short_string_unchanged(self):
        assert softwrap("ATCG", 80) == "ATCG"

    def test_wraps_at_specified_width(self):
        result = softwrap("A" * 20, 10)
        lines = result.split("\n")
        assert len(lines) == 2
        assert all(len(line) == 10 for line in lines)

    def test_exact_width_no_wrap(self):
        assert softwrap("A" * 10, 10) == "A" * 10

    def test_one_over_width_wraps(self):
        result = softwrap("A" * 11, 10)
        assert "\n" in result
        assert result.split("\n")[0] == "A" * 10

    def test_empty_string(self):
        assert softwrap("", 80) == ""

    def test_default_width_is_60(self):
        assert softwrap("A" * 61).split("\n") == ["A" * 60, "A"]


class TestWriteFasta:
    def test_wraps_at_60(self):
        buf = io.StringIO()
        write_fasta(buf, "ctg1 desc", "A" * 61)
        assert buf.getvalue() == ">ctg1 desc\n" + "A" * 60 + "\nA\n"

    def test_empty_sequence_has_no_blank_line(self):
        buf = io.StringIO()
        write_fasta(buf, "empty", "")
        assert buf.getvalue() == ">empty\n"


# ---------------------------------------------------------------------------
# check_file
# ---------------------------------------------------------------------------


class TestCheckfile:
    def test_existing_nonempty_file_is_true(self, fasta_file):
        assert check_file(str(fasta_file)) is True

    def test_empty_file_is_false(self, empty_file):
        assert check_file(str(empty_file)) is False

    def test_missing_file_is_false(self, tmp_path):
        assert check_file(str(tmp_path / "nonexistent.fa")) is False

    def test_symlink_to_nonempty_file_is_true(self, tmp_path, fasta_file):
        link = tmp_path / "link.fa"
        link.symlink_to(fasta_file)
        assert check_file(str(link)) is True

    def test_broken_symlink_is_false(self, tmp_path):
        link = tmp_path / "broken_link.fa"
        link.symlink_to(tmp_path / "nonexistent_target.fa")
        assert check_file(str(link)) is False


# ---------------------------------------------------------------------------
# fasta_stats / filter_fasta
# ---------------------------------------------------------------------------


class TestFastastats:
    def test_count(self, fasta_file):
        count, total_len = fasta_stats(str(fasta_file))
        assert count == 3

    def test_total_length(self, fasta_file):
        # SEQ1=20, SEQ2=80, SEQ3=150 → 250
        count, total_len = fasta_stats(str(fasta_file))
        assert total_len == 250

    def test_wrapped_sequence_length(self, tmp_path):
        p = tmp_path / "wrapped.fa"
        p.write_text(">a desc\nACGT\nAC\n>b\nGG\n")
        assert fasta_stats(str(p)) == (2, 8)


class TestFilterFasta:
    def test_drops_by_id_and_returns_stats(self, tmp_path):
        src = tmp_path / "in.fa"
        src.write_text(">a desc\nACGT\nAC\n>b\nGG\n>c\nTTT\n")
        out = tmp_path / "out.fa"
        assert filter_fasta(str(src), str(out), lambda seq_id: seq_id != "b") == (2, 9)
        assert out.read_text() == ">a desc\nACGTAC\n>c\nTTT\n"

    def test_matches_biopython_seqio_write(self, tmp_path):
        from Bio import SeqIO

        src = tmp_path / "in.fa"
        src.write_text(">a desc\n" + "ACGT" * 40 + "\n>b\n" + "G" * 60 + "\n>c\n" + "T" * 61 + "\n")
        out = tmp_path / "out.fa"
        expected = tmp_path / "expected.fa"
        filter_fasta(str(src), str(out), lambda seq_id: True)
        with open(src) as fh:
            SeqIO.write(SeqIO.parse(fh, "fasta"), str(expected), "fasta")
        assert out.read_text() == expected.read_text()

    def test_keep_none(self, tmp_path):
        src = tmp_path / "in.fa"
        src.write_text(">a\nACGT\n")
        out = tmp_path / "out.fa"
        assert filter_fasta(str(src), str(out), lambda seq_id: False) == (0, 0)
        assert out.read_text() == ""


# ---------------------------------------------------------------------------
# count_fastq
# ---------------------------------------------------------------------------


class TestCountFastq:
    def test_plain_fastq(self, fastq_file):
        assert count_fastq(str(fastq_file)) == 10

    def test_gzip_fastq(self, gz_fastq_file):
        assert count_fastq(str(gz_fastq_file)) == 10

    def test_single_read(self, tmp_path):
        p = tmp_path / "one.fastq"
        p.write_text(make_fastq_text(1))
        assert count_fastq(str(p)) == 1

    def test_gz_50_reads(self, tmp_path):
        p = tmp_path / "fifty.fastq.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(make_fastq_text(50))
        assert count_fastq(str(p)) == 50

    def test_missing_trailing_newline(self, tmp_path):
        p = tmp_path / "nonl.fastq"
        p.write_text(make_fastq_text(3).rstrip("\n"))
        assert count_fastq(str(p)) == 3

    def test_gz_without_pigz_or_gzip(self, gz_fastq_file, monkeypatch):
        monkeypatch.setattr("aaftf.utility.shutil.which", lambda name: None)
        assert count_fastq(str(gz_fastq_file)) == 10

    def test_bad_gzip_raises(self, tmp_path):
        p = tmp_path / "bad.fastq.gz"
        p.write_bytes(b"not a valid gzip file")
        with pytest.raises((OSError, subprocess.CalledProcessError)):
            count_fastq(str(p))


# ---------------------------------------------------------------------------
# safe_remove
# ---------------------------------------------------------------------------


class TestSafeRemove:
    def test_remove_file(self, tmp_path):
        p = tmp_path / "todelete.txt"
        p.write_text("bye")
        safe_remove(str(p))
        assert not p.exists()

    def test_remove_directory(self, tmp_path):
        d = tmp_path / "subdir"
        d.mkdir()
        (d / "file.txt").write_text("hi")
        safe_remove(str(d))
        assert not d.exists()

    def test_remove_nonexistent_is_noop(self, tmp_path):
        # Should not raise
        safe_remove(str(tmp_path / "ghost"))

    def test_symlink_to_directory_removes_link_not_target(self, tmp_path):
        target = tmp_path / "target"
        target.mkdir()
        (target / "keep.txt").write_text("keep")
        link = tmp_path / "link"
        link.symlink_to(target)
        safe_remove(str(link))
        assert not link.exists() and not link.is_symlink()
        assert (target / "keep.txt").exists()

    def test_broken_symlink_removed(self, tmp_path):
        link = tmp_path / "broken"
        link.symlink_to(tmp_path / "missing")
        safe_remove(str(link))
        assert not link.is_symlink()


class TestCleanupWorkdir:
    def _workdir(self, tmp_path):
        d = tmp_path / "wd"
        d.mkdir()
        return d

    def test_removes_auto_generated_workdir(self, tmp_path):
        d = self._workdir(tmp_path)
        cleanup_workdir(str(d), debug=False, custom_workdir=False)
        assert not d.exists()

    def test_keeps_user_supplied_workdir(self, tmp_path):
        d = self._workdir(tmp_path)
        cleanup_workdir(str(d), debug=False, custom_workdir=True)
        assert d.exists()

    def test_keeps_workdir_when_debugging(self, tmp_path):
        d = self._workdir(tmp_path)
        cleanup_workdir(str(d), debug=True, custom_workdir=False)
        assert d.exists()


# ---------------------------------------------------------------------------
# samtools helpers
# ---------------------------------------------------------------------------


class TestSamtoolsSortCmd:
    def test_basic(self):
        assert samtools_sort_cmd("-", "out.bam", 4) == ["samtools", "sort", "-@", "4", "-o", "out.bam", "-"]

    def test_memory_and_tmp_prefix(self):
        cmd = samtools_sort_cmd("in.sam", "out.bam", 2, memory_per_thread="1G", tmp_prefix="tmp")
        assert cmd == ["samtools", "sort", "-@", "2", "-m", "1G", "-o", "out.bam", "-T", "tmp", "in.sam"]

    def test_write_index_requests_bai(self):
        cmd = samtools_sort_cmd("-", "out.bam", 1, write_index=True)
        assert cmd == ["samtools", "sort", "-@", "1", "--write-index", "-o", "out.bam##idx##out.bam.bai", "-"]


_SAM = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\nr1\t0\tc1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
# 2 reads: r1 mapped (plus a supplementary and a secondary record), r2 unmapped
_SAM_MIXED = _SAM + "r1\t2048\tc1\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\nr1\t256\tc1\t30\t0\t4M\t*\t0\t0\tACGT\tIIII\nr2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")
class TestAlignToSortedBam:
    def test_writes_bam_relative_to_current_dir_not_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sub = tmp_path / "work"
        sub.mkdir()
        (sub / "in.sam").write_text(_SAM)
        align_to_sorted_bam(["cat", "in.sam"], "out.bam", cwd=str(sub))
        assert (tmp_path / "out.bam").stat().st_size > 0
        assert (tmp_path / "out.bam.bai").exists()
        assert not (sub / "out.bam").exists()

    def test_bam_read_count_ignores_secondary_and_supplementary(self, tmp_path):
        (tmp_path / "in.sam").write_text(_SAM_MIXED)
        bam = tmp_path / "out.bam"
        align_to_sorted_bam(["cat", str(tmp_path / "in.sam")], str(bam))
        assert bam_read_count(str(bam)) == (1, 1)

    def test_failing_aligner_raises_and_removes_partial_bam(self, tmp_path):
        bam = tmp_path / "out.bam"
        with pytest.raises(RuntimeError):
            align_to_sorted_bam(["false"], str(bam))
        assert not bam.exists()


class TestRunCmd:
    def test_prints_and_returns_result(self, capsys):
        result = run_cmd(["true"])
        assert result.returncode == 0
        assert "true" in capsys.readouterr().out

    def test_stderr_hidden_unless_debug(self, capfd):
        run_cmd(["sh", "-c", "echo oops >&2"])
        assert "oops" not in capfd.readouterr().err
        run_cmd(["sh", "-c", "echo oops >&2"], debug=True)
        assert "oops" in capfd.readouterr().err

    def test_quiet_stdout(self, capfd):
        run_cmd(["echo", "hello"], quiet_stdout=True)
        # only the printed "CMD: echo hello" line should appear, not echo's own output
        assert "hello" not in capfd.readouterr().out.splitlines()

    def test_stdout_to_file(self, tmp_path):
        out = tmp_path / "out.txt"
        with open(out, "w") as fh:
            run_cmd(["echo", "hello"], stdout=fh)
        assert out.read_text() == "hello\n"


class TestExecute:
    def test_yields_lines_and_prints_command(self, capsys):
        assert list(execute(["printf", "a\\nb\\n"])) == ["a\n", "b\n"]
        assert "printf" in capsys.readouterr().out

    def test_quiet_does_not_print(self, capsys):
        list(execute(["true"], quiet=True))
        assert capsys.readouterr().out == ""

    def test_cwd(self, tmp_path):
        assert list(execute(["pwd"], cwd=str(tmp_path), quiet=True)) == [f"{tmp_path}\n"]

    def test_failure_raises_after_output(self):
        with pytest.raises(subprocess.CalledProcessError):
            list(execute(["sh", "-c", "echo x; exit 3"], quiet=True))

    def test_stopping_early_ends_the_command_without_raising(self):
        gen = execute(["yes"], quiet=True)
        assert next(gen) == "y\n"
        gen.close()  # what a `break` in the caller's loop does


class TestNextStepName:
    def test_splits_at_first_underscore(self):
        assert next_step_name("strain_vecscreen.fasta", ".taxpurge.fasta") == "strain.taxpurge.fasta"

    def test_splits_at_first_dot_without_underscore(self):
        assert next_step_name("strain.rmdup.fasta", ".polish.fasta") == "strain.polish.fasta"

    def test_no_separator(self):
        assert next_step_name("strain", ".final.fasta") == "strain.final.fasta"


class TestBasenameFromReads:
    def test_underscore(self):
        assert basename_from_reads("/data/Sample1_R1.fastq.gz") == "Sample1"

    def test_dot(self):
        assert basename_from_reads("Sample1.fq.gz") == "Sample1"

    def test_uses_file_name_not_directory(self):
        assert basename_from_reads("/my_dir/reads") == "reads"


class TestMakeWorkdir:
    def test_auto_named_and_created(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        workdir, custom = make_workdir(None, "filter")
        assert custom is False
        assert workdir.startswith("aaftf-filter_") and (tmp_path / workdir).is_dir()

    def test_user_supplied_is_created_and_flagged(self, tmp_path):
        target = tmp_path / "a" / "b"
        workdir, custom = make_workdir(str(target), "filter")
        assert custom is True and workdir == str(target) and target.is_dir()


class TestRequireTools:
    def test_present_tools_pass(self):
        require_tools(["sh"])

    def test_missing_tool_raises(self):
        with pytest.raises(FileNotFoundError) as exc:
            require_tools(["sh", "definitely_not_a_tool_xyz"])
        assert "definitely_not_a_tool_xyz" in str(exc.value)


class TestOpenAndConcat:
    def test_open_maybe_gz_reads_both(self, tmp_path):
        plain, gz = tmp_path / "a.txt", tmp_path / "b.txt.gz"
        plain.write_text("plain\n")
        with gzip.open(gz, "wt") as fh:
            fh.write("zipped\n")
        with open_maybe_gz(plain) as fh:
            assert fh.read() == "plain\n"
        with open_maybe_gz(gz) as fh:
            assert fh.read() == "zipped\n"

    def test_concat_files_mixes_plain_and_gz(self, tmp_path):
        plain, gz, out = tmp_path / "a.fa", tmp_path / "b.fa.gz", tmp_path / "all.fa"
        plain.write_text(">a\nAC\n")
        with gzip.open(gz, "wt") as fh:
            fh.write(">b\nGT\n")
        concat_files([str(plain), str(gz)], str(out))
        assert out.read_text() == ">a\nAC\n>b\nGT\n"


class TestDatabaseFolders:
    def test_unset_uses_home_cache(self, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
        assert home_db_cache() == (tmp_path / "aaftf").resolve()
        assert db_dirs() == [home_db_cache()]

    def test_home_cache_without_xdg(self, monkeypatch, tmp_path):
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        monkeypatch.setenv("HOME", str(tmp_path))
        assert home_db_cache() == (tmp_path / ".cache" / "aaftf").resolve()

    def test_path_like_list(self, monkeypatch, tmp_path):
        shared, mine = tmp_path / "shared", tmp_path / "mine"
        monkeypatch.setenv("AAFTF_DB", f"{shared}{os.pathsep}{mine}")
        assert db_dirs() == [shared.resolve(), mine.resolve()]

    def test_finds_file_in_any_folder_first_wins(self, monkeypatch, tmp_path):
        shared, mine = tmp_path / "shared", tmp_path / "mine"
        shared.mkdir()
        mine.mkdir()
        (mine / "UniVec").write_text("mine")
        monkeypatch.setenv("AAFTF_DB", f"{shared}{os.pathsep}{mine}")
        assert find_db_file("UniVec") == str(mine.resolve() / "UniVec")
        (shared / "UniVec").write_text("shared")
        assert find_db_file("UniVec") == str(shared.resolve() / "UniVec")
        assert find_db_file("missing") is None

    def test_find_db_index_needs_a_complete_rocksdb_folder(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        index = tmp_path / "db.rocksdb"
        index.mkdir()
        assert find_db_index("db.rocksdb") is None  # no CURRENT file: empty or partly copied
        (index / "CURRENT").write_text("MANIFEST-000005\n")
        assert find_db_index("db.rocksdb") == str(tmp_path.resolve() / "db.rocksdb")
        assert find_db_index("missing.rocksdb") is None

    def test_matches_with_or_without_gz(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        (tmp_path / "contam_in_euks.fa").write_text(">a\n")
        (tmp_path / "UniVec.gz").write_bytes(b"")
        assert find_db_file("contam_in_euks.fa.gz") == str(tmp_path.resolve() / "contam_in_euks.fa")
        assert find_db_file("UniVec") == str(tmp_path.resolve() / "UniVec.gz")
        (tmp_path / "UniVec").write_text(">b\n")
        assert find_db_file("UniVec") == str(tmp_path.resolve() / "UniVec")

    def test_missing_file_goes_to_first_writable_folder(self, monkeypatch, tmp_path):
        readonly, mine = tmp_path / "readonly", tmp_path / "mine"
        readonly.mkdir()
        readonly.chmod(0o555)
        try:
            monkeypatch.setenv("AAFTF_DB", f"{readonly}{os.pathsep}{mine}")
            assert db_file("UniVec") == str(mine.resolve() / "UniVec")
            assert mine.is_dir()
        finally:
            readonly.chmod(0o755)

    def test_force_ignores_existing_copy(self, monkeypatch, tmp_path):
        shared, mine = tmp_path / "shared", tmp_path / "mine"
        shared.mkdir()
        (shared / "UniVec").write_text("shared")
        shared.chmod(0o555)
        try:
            monkeypatch.setenv("AAFTF_DB", f"{shared}{os.pathsep}{mine}")
            assert db_file("UniVec") == str(shared.resolve() / "UniVec")
            assert db_file("UniVec", force=True) == str(mine.resolve() / "UniVec")
        finally:
            shared.chmod(0o755)

    def test_no_writable_folder_falls_back_to_home_cache(self, monkeypatch, tmp_path):
        readonly = tmp_path / "readonly"
        readonly.mkdir()
        readonly.chmod(0o555)
        try:
            monkeypatch.setenv("AAFTF_DB", str(readonly))
            assert db_write_dir() == home_db_cache()
        finally:
            readonly.chmod(0o755)


class TestHomeCacheWarning:
    def test_warns_once_when_unset(self, caplog):
        warn_if_home_cache()
        warn_if_home_cache()
        assert caplog.text.count("AAFTF_DB is not set") == 1
        assert "export AAFTF_DB=" in caplog.text

    def test_no_warning_when_aaftf_db_writable(self, monkeypatch, tmp_path, caplog):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        warn_if_home_cache()
        assert "AAFTF_DB" not in caplog.text

    def test_download_into_home_cache_warns(self, tmp_path, caplog):
        src = tmp_path / "src.txt"
        src.write_text("data")
        download_file(src.as_uri(), db_file("src.txt"))
        assert "AAFTF_DB is not set" in caplog.text


class TestDownloadFile:
    def test_downloads_and_returns_dest(self, tmp_path):
        src = tmp_path / "src.txt"
        src.write_text("data")
        dest = tmp_path / "sub" / "dest.txt"
        assert download_file(src.as_uri(), str(dest)) == str(dest)
        assert dest.read_text() == "data"

    def test_existing_file_not_redownloaded(self, tmp_path):
        dest = tmp_path / "dest.txt"
        dest.write_text("cached")
        download_file((tmp_path / "missing.txt").as_uri(), str(dest))
        assert dest.read_text() == "cached"

    def test_failed_download_leaves_no_partial_file(self, tmp_path):
        dest = tmp_path / "dest.txt"
        with pytest.raises(Exception):
            download_file((tmp_path / "missing.txt").as_uri(), str(dest))
        assert not dest.exists() and not (tmp_path / "dest.txt.tmp").exists()


class TestSetupLogging:
    def _emit(self, capsys, **kwargs):
        import logging

        setup_logging(**kwargs)
        log = logging.getLogger("aaftf.test")
        log.debug("dbg")
        log.info("inf")
        log.warning("wrn")
        return capsys.readouterr().err

    def test_default_shows_info_and_prefixes_warnings(self, capsys):
        err = self._emit(capsys)
        assert "inf" in err and "dbg" not in err
        assert "WARNING: wrn" in err

    def test_quiet_shows_only_warnings(self, capsys):
        err = self._emit(capsys, quiet=True)
        assert "inf" not in err and "WARNING: wrn" in err

    def test_debug_shows_debug(self, capsys):
        assert "dbg" in self._emit(capsys, debug=True)

    def test_no_color_when_not_a_terminal(self, capsys):
        assert "\033[" not in self._emit(capsys)


class TestAvailableCpus:
    def test_slurm_allocation_wins(self, monkeypatch):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "3")
        assert available_cpus() == 3

    @pytest.mark.parametrize("value", ["0", "four", "-2"])
    def test_invalid_slurm_value_ignored(self, monkeypatch, value):
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", value)
        monkeypatch.setattr("aaftf.utility.os.sched_getaffinity", lambda pid: {0, 1, 2}, raising=False)
        assert available_cpus() == 3


class TestGetRam:
    """get_ram() rounds available RAM down (to a multiple of 10 above 10 GB, else 1 decimal minus 0.1)."""

    @staticmethod
    def _ram(monkeypatch, avail_gb, max_lim=0):
        from types import SimpleNamespace

        from aaftf.utility import get_ram

        monkeypatch.setattr("aaftf.utility.psutil.virtual_memory", lambda: SimpleNamespace(available=avail_gb * 1024**3))
        return get_ram(max_lim)

    @pytest.mark.parametrize("avail, expected", [(25.7, 20.0), (10.5, 10.0), (3.05, 2.9), (0.04, 0.0)])
    def test_rounds_down(self, monkeypatch, avail, expected):
        result = self._ram(monkeypatch, avail)
        assert result == pytest.approx(expected)
        assert result <= avail

    def test_max_lim_caps(self, monkeypatch):
        assert self._ram(monkeypatch, 25.7, max_lim=8) == 8

    def test_max_lim_above_available_does_not_raise_it(self, monkeypatch):
        assert self._ram(monkeypatch, 25.7, max_lim=64) == 20.0


class TestDbDirs:
    def test_empty_entries_fall_back_to_home_cache(self, monkeypatch):
        from aaftf.utility import db_dirs, home_db_cache

        monkeypatch.setenv("AAFTF_DB", os.pathsep)
        assert db_dirs() == [home_db_cache()]

    def test_tilde_expanded(self, monkeypatch, tmp_path):
        from aaftf.utility import db_dirs

        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("AAFTF_DB", f"~/dbs{os.pathsep}{os.pathsep}")
        assert db_dirs() == [(tmp_path / "dbs").resolve()]

    def test_write_dir_skips_folder_that_cannot_be_created(self, monkeypatch, tmp_path):
        from pathlib import Path

        from aaftf.utility import db_write_dir

        bad, good = (tmp_path / "bad").resolve(), (tmp_path / "good").resolve()
        monkeypatch.setenv("AAFTF_DB", f"{bad}{os.pathsep}{good}")
        real_mkdir = Path.mkdir

        def _mkdir(self, *args, **kwargs):
            if self == bad:
                raise PermissionError(str(self))
            return real_mkdir(self, *args, **kwargs)

        monkeypatch.setattr(Path, "mkdir", _mkdir)
        assert db_write_dir() == good
        assert good.is_dir() and not bad.exists()


def test_make_workdir_without_setup_logging_writes_no_log(tmp_path):
    from aaftf.utility import cleanup_workdir, make_workdir

    folder, custom = make_workdir(str(tmp_path / "wd"), "step")
    assert custom and (tmp_path / "wd").is_dir()
    cleanup_workdir(folder, False, custom)
    assert not (tmp_path / "wd" / "step.log").exists()

    def test_uses_affinity_without_slurm(self, monkeypatch):
        monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
        monkeypatch.setattr("aaftf.utility.os.sched_getaffinity", lambda pid: {0, 1}, raising=False)
        assert available_cpus() == 2

    def test_falls_back_to_cpu_count_without_affinity(self, monkeypatch):
        monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
        monkeypatch.delattr("aaftf.utility.os.sched_getaffinity", raising=False)
        monkeypatch.setattr("aaftf.utility.os.cpu_count", lambda: 5)
        assert available_cpus() == 5


class TestRequireDatabases:
    def test_returns_paths_of_stored_databases(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        (tmp_path / "UniVec").write_text(">v\nACGT\n")
        assert require_databases(["univec"]) == [str(tmp_path.resolve() / "UniVec")]

    def test_missing_raises_with_download_command(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        (tmp_path / "UniVec").write_text(">v\nACGT\n")
        with pytest.raises(FileNotFoundError) as exc:
            require_databases(["univec", "euks", "proks"], hint="or pass --x")
        assert "missing database(s): euks, proks" in str(exc.value)
        assert "AAFTF database euks proks (or pass --x)" in str(exc.value)

    def test_never_downloads(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        with patch("aaftf.utility.download_file") as download:
            with pytest.raises(FileNotFoundError):
                require_databases(["univec"])
        download.assert_not_called()


class TestPafHits:
    def test_parses_named_fields_and_skips_short_lines(self):
        paf = "q1\t1000\t10\t990\t-\tt1\t5000\t100\t1080\t950\t980\t60\ttp:A:P\nnot paf\n"
        hits = list(paf_hits(["printf", paf], quiet=True))
        assert hits == [PafHit("q1", 1000, 10, 990, "-", "t1", 5000, 100, 1080, 950, 980, 60)]
        assert hits[0].strand == "-" and hits[0].target_end == 1080


def _rev_comp(seq):
    """Reverse complement the way callers use COMPLEMENT."""
    return seq.translate(COMPLEMENT)[::-1]


class TestComplement:
    """COMPLEMENT (with slicing) reverse complements DNA, keeping case and IUPAC codes."""

    def test_simple(self):
        assert _rev_comp("ATCG") == "CGAT"

    def test_complement_only(self):
        assert _rev_comp("AAAA") == "TTTT"
        assert _rev_comp("CCCC") == "GGGG"

    def test_palindrome(self):
        assert _rev_comp("AATTAATT") == "AATTAATT"

    def test_preserves_case(self):
        assert _rev_comp("atcg") == "cgat"
        assert _rev_comp("AAcgTT") == "AAcgTT"

    def test_iupac_codes(self):
        assert _rev_comp("RYKMN") == "NKMRY"

    def test_longer_sequence(self):
        assert _rev_comp("ATCGATCG") == "CGATCGAT"

    def test_all_bases(self):
        # A<->T, C<->G
        assert _rev_comp("ACGT") == "ACGT"
