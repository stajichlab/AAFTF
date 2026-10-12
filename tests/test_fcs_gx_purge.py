"""Unit tests for aaftf/fcs_gx_purge.py (run_gx.py, the filesystem type and the speed test are mocked)."""

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf import fcs_gx_purge
from aaftf.fcs_gx_purge import cleaned_name, read_latency_ms, run

pytestmark = pytest.mark.unit

_HEADER = "##[header]\n#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n"


@pytest.fixture
def scratch(tmp_path, monkeypatch):
    """Return a --copy_to folder; report plenty of RAM, a network filesystem and slow reads (so ``auto`` copies)."""
    monkeypatch.setattr(fcs_gx_purge, "get_ram", lambda max_lim=0: 1000.0)
    monkeypatch.setattr(fcs_gx_purge, "_fs_type", lambda path: "lustre")
    monkeypatch.setattr(fcs_gx_purge, "read_latency_ms", lambda path: 5.0)
    monkeypatch.setattr(fcs_gx_purge, "_same_filesystem", lambda path, folder: False)
    monkeypatch.setattr("aaftf.fcs_gx_purge.shutil.which", lambda name: f"/usr/bin/{name}")
    return tmp_path / "scratch"


@pytest.fixture
def db(tmp_path):
    """A fake FCS-GX database prefix with a .gxi and a .gxs file."""
    db = tmp_path / "gxdb" / "all"
    db.parent.mkdir()
    Path(f"{db}.gxi").write_text("index")
    Path(f"{db}.gxs").write_text("sequences")
    return db


def _fake_run_gx(rows=(), returncode=0, write_report=True, seen=None):
    """Return a subprocess.run stand-in that writes an FCS-GX report with ``rows`` into --out-dir."""

    def _run(cmd, **kwargs):
        if seen is not None:
            gx_db = cmd[cmd.index("--gx-db") + 1]
            seen.append((cmd[cmd.index("--fasta") + 1], gx_db, sorted(p.name for p in Path(gx_db).parent.iterdir())))
        if write_report:
            stem = Path(cmd[cmd.index("--fasta") + 1]).stem
            Path(cmd[cmd.index("--out-dir") + 1], f"{stem}.4890.fcs_gx_report.txt").write_text(_HEADER + "".join("\t".join(r) + "\n" for r in rows))
        return subprocess.CompletedProcess(cmd, returncode)

    return _run


def _asm(path, ids=("keep", "bad", "trimme")):
    path.write_text("".join(f">{i}\nACGT\n" for i in ids))
    return path


def _ids(path):
    return [line[1:] for line in Path(path).read_text().split() if line.startswith(">")]


class TestCleanedName:
    @pytest.mark.parametrize("name", ["s.final.fasta", "s.final.fa", "s.final.fna.gz", "s.final"])
    def test_suffix_before_fasta(self, name):
        assert cleaned_name(f"/a/{name}") == "/a/s.final.fcsgx_cleaned.fasta"

    def test_outdir(self):
        assert cleaned_name("/a/s.fasta", "/out") == "/out/s.fcsgx_cleaned.fasta"


class TestScreen:
    def test_only_exclude_rows_are_dropped(self, tmp_path, scratch, db, caplog):
        asm = _asm(tmp_path / "asm.fasta")
        rows = [["bad", "1", "4", "4", "EXCLUDE", "anml:insects", "100", "fly"], ["trimme", "1", "2", "4", "TRIM", "prok:bact", "50", "E. coli"]]
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx(rows)):
            run(taxid=4890, input=[str(asm)], db=str(db), workdir=str(tmp_path / "wd"), copy_to=str(scratch))
        assert _ids(tmp_path / "asm.fcsgx_cleaned.fasta") == ["keep", "trimme"]
        assert "TRIM for trimme:1..2" in caplog.text
        assert "EXCLUDE" in (tmp_path / "asm.fcs_gx_report.txt").read_text()

    def test_action_prefix_accepted(self, tmp_path, scratch, db):
        asm = _asm(tmp_path / "asm.fasta")
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx([["bad", "1", "4", "4", "ACTION_EXCLUDE", "x", "100", "y"]])):
            run(taxid=4890, input=[str(asm)], db=str(db), copy_to=str(scratch))
        assert _ids(tmp_path / "asm.fcsgx_cleaned.fasta") == ["keep", "trimme"]

    def test_failed_run_gx_raises_and_removes_copy(self, tmp_path, scratch, db):
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx(returncode=1)), pytest.raises(RuntimeError, match="run_gx.py failed"):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch))
        assert not list(scratch.iterdir())

    def test_missing_report_raises(self, tmp_path, scratch, db):
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx(write_report=False)), pytest.raises(RuntimeError, match="did not write"):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch))


class TestCopy:
    """Where run_gx.py reads the database from: in place, or one per-job copy under --copy_to."""

    def _db_used(self, tmp_path, db, **opts):
        seen = []
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx(seen=seen)):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), **opts)
        return seen[0]

    def test_slow_db_copied_once_for_all_assemblies(self, tmp_path, scratch, db, monkeypatch):
        monkeypatch.setenv("SLURM_JOB_ID", "777")
        a, b = _asm(tmp_path / "a.fasta"), _asm(tmp_path / "b.fa", ids=("x", "bad"))
        seen = []
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx([["bad", "1", "4", "4", "EXCLUDE", "x", "100", "y"]], seen=seen)):
            run(taxid=4890, input=[str(a), str(b)], db=str(db), outdir=str(tmp_path / "out"), copy_to=str(scratch))
        assert _ids(tmp_path / "out/a.fcsgx_cleaned.fasta") == ["keep", "trimme"]
        assert _ids(tmp_path / "out/b.fcsgx_cleaned.fasta") == ["x"]
        (_, db_a, files), (_, db_b, _) = seen
        assert db_a == db_b and Path(db_a).parent.parent == scratch and Path(db_a).name == "all"
        assert Path(db_a).parent.name.startswith("aaftf-gxdb-777-")
        assert files == ["all.gxi", "all.gxs"]
        assert not list(scratch.iterdir())  # copy removed afterwards

    def test_tmpfs_db_used_in_place_without_speed_test_or_ram_check(self, tmp_path, scratch, db, monkeypatch):
        monkeypatch.setattr(fcs_gx_purge, "_fs_type", lambda path: "tmpfs")
        monkeypatch.setattr(fcs_gx_purge, "get_ram", lambda max_lim=0: 0.0)
        monkeypatch.setattr(fcs_gx_purge, "read_latency_ms", lambda path: pytest.fail("speed test run"))
        assert self._db_used(tmp_path, db, copy_to=str(scratch))[1] == str(db)
        assert not scratch.exists()

    def test_fast_db_used_in_place(self, tmp_path, scratch, db, monkeypatch):
        monkeypatch.setattr(fcs_gx_purge, "read_latency_ms", lambda path: 0.05)
        assert self._db_used(tmp_path, db, copy_to=str(scratch))[1] == str(db)
        assert Path(f"{db}.gxi").is_file()

    def test_speed_test_reads_largest_file(self, tmp_path, scratch, db, monkeypatch):
        tested = []
        monkeypatch.setattr(fcs_gx_purge, "read_latency_ms", lambda path: tested.append(path) or 0.05)
        self._db_used(tmp_path, db, copy_to=str(scratch))
        assert tested == [f"{db}.gxs"]

    def test_copy_defaults_to_tmpdir(self, tmp_path, scratch, db, monkeypatch):
        monkeypatch.setenv("TMPDIR", str(scratch))
        assert Path(self._db_used(tmp_path, db)[1]).parent.parent == scratch

    def test_copy_falls_back_without_tmpdir(self, tmp_path, scratch, db, monkeypatch):
        monkeypatch.delenv("TMPDIR", raising=False)
        monkeypatch.setattr(fcs_gx_purge, "FALLBACK_COPY_TO", str(scratch))
        assert Path(self._db_used(tmp_path, db)[1]).parent.parent == scratch


class TestDbArgument:
    """-d/--db takes the database path prefix, its .gxi file, or the folder holding it."""

    @pytest.mark.parametrize("form", ["prefix", "gxi", "folder"])
    def test_forms(self, tmp_path, scratch, db, monkeypatch, form):
        monkeypatch.setattr(fcs_gx_purge, "read_latency_ms", lambda path: 0.05)
        arg = {"prefix": str(db), "gxi": f"{db}.gxi", "folder": str(db.parent)}[form]
        seen = []
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx(seen=seen)):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=arg, copy_to=str(scratch))
        assert seen[0][1] == str(db)

    def test_folder_with_two_indexes(self, tmp_path, scratch, db):
        Path(db.parent, "other.gxi").write_text("")
        with pytest.raises(FileNotFoundError, match="exactly one"):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db.parent), copy_to=str(scratch))

    def test_all_prefix_files_copied(self, tmp_path, scratch, db):
        Path(f"{db}.taxa.tsv").write_text("t")
        Path(db.parent, "unrelated.txt").write_text("u")
        seen = []
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx(seen=seen)):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch))
        assert seen[0][2] == ["all.gxi", "all.gxs", "all.taxa.tsv"]


class TestGuards:
    """Checks that run before the database is copied or run_gx.py is started."""

    def test_not_enough_ram(self, tmp_path, scratch, db, monkeypatch):
        monkeypatch.setattr(fcs_gx_purge, "get_ram", lambda max_lim=0: 0.0)
        with patch("aaftf.fcs_gx_purge.subprocess.run") as gx, pytest.raises(MemoryError, match="sbatch --mem="):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch))
        gx.assert_not_called()
        assert not scratch.exists()

    def test_not_enough_copy_space(self, tmp_path, scratch, db, monkeypatch):
        from types import SimpleNamespace

        monkeypatch.setattr("aaftf.fcs_gx_purge.shutil.disk_usage", lambda path: SimpleNamespace(free=0))
        with patch("aaftf.fcs_gx_purge.subprocess.run") as gx, pytest.raises(OSError, match="--copy_to"):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch))
        gx.assert_not_called()

    def test_missing_db(self, tmp_path, scratch):
        with pytest.raises(FileNotFoundError, match="all.gxi"):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(tmp_path / "nodb" / "all"), copy_to=str(scratch))

    def test_missing_input(self, tmp_path, scratch, db):
        with pytest.raises(FileNotFoundError, match="missing.fasta"):
            run(taxid=4890, input=[str(tmp_path / "missing.fasta")], db=str(db), copy_to=str(scratch))

    def test_clashing_outputs(self, tmp_path, scratch, db):
        (tmp_path / "x").mkdir()
        a, b = _asm(tmp_path / "s.fasta"), _asm(tmp_path / "x" / "s.fa")
        with pytest.raises(ValueError, match="same"):
            run(taxid=4890, input=[str(a), str(b)], db=str(db), outdir=str(tmp_path / "out"), copy_to=str(scratch))


class TestReadLatency:
    def test_measures_a_real_file(self, tmp_path):
        f = tmp_path / "blob"
        f.write_bytes(b"\0" * 1_000_000)
        assert read_latency_ms(str(f), reads=50) > 0

    def test_stops_at_time_limit(self, tmp_path):
        f = tmp_path / "blob"
        f.write_bytes(b"\0" * 10_000)
        with patch("aaftf.fcs_gx_purge.os.pread") as pread:
            read_latency_ms(str(f), reads=10_000, max_seconds=0)
        assert pread.call_count == 1


class TestFsType:
    def test_longest_mount_wins(self, tmp_path):
        mounts = "/dev/sda1 / xfs rw 0 0\ntmpfs /dev/shm tmpfs rw 0 0\nfs /scratch lustre rw 0 0\n"
        with patch("builtins.open", lambda *a, **k: __import__("io").StringIO(mounts)), patch("aaftf.fcs_gx_purge.os.path.realpath", lambda p: p):
            from aaftf.fcs_gx_purge import _fs_type

            assert _fs_type("/dev/shm/gxdb/all") == "tmpfs"
            assert _fs_type("/scratch/gxdb/all") == "lustre"
            assert _fs_type("/home/u/gxdb") == "xfs"
            assert _fs_type("/scratchy/x") == "xfs"


def test_cli_options(tmp_path, scratch, db):
    import sys

    from aaftf.main import main

    a, b = _asm(tmp_path / "a.fasta"), _asm(tmp_path / "b.fasta")
    argv = ["AAFTF", "fcs_gx_purge", "-i", str(a), str(b), "-d", str(db), "-t", "4890", "-c", "3", "--copy_to", str(scratch)]
    with patch.object(sys, "argv", argv), patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx()):
        assert main() == 0
    assert (tmp_path / "a.fcsgx_cleaned.fasta").is_file() and (tmp_path / "b.fcsgx_cleaned.fasta").is_file()


def test_cli_requires_taxid(tmp_path, db):
    import sys

    from aaftf.main import main

    with patch.object(sys, "argv", ["AAFTF", "fcs_gx_purge", "-i", "a.fasta", "-d", str(db)]), pytest.raises(SystemExit):
        main()


class TestRunGxSetup:
    def test_missing_run_gx_fails_before_copy(self, tmp_path, scratch, db, monkeypatch):
        monkeypatch.setattr("aaftf.fcs_gx_purge.shutil.which", lambda name: None)
        with pytest.raises(FileNotFoundError, match="run_gx.py"):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch))
        assert not scratch.exists()

    @pytest.mark.parametrize("cpus, expected", [(5, "5"), (None, "7")])
    def test_cpus_passed_as_gx_num_cores(self, tmp_path, scratch, db, monkeypatch, cpus, expected):
        monkeypatch.setattr(fcs_gx_purge, "available_cpus", lambda: 7)
        envs = []
        fake = _fake_run_gx()

        def _run(cmd, **kwargs):
            envs.append(kwargs["env"])
            return fake(cmd, **kwargs)

        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_run):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch), cpus=cpus)
        assert envs[0]["GX_NUM_CORES"] == expected


class TestCopyTarget:
    def _db_used(self, tmp_path, db, scratch):
        seen = []
        with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run_gx(seen=seen)):
            run(taxid=4890, input=[str(_asm(tmp_path / "asm.fasta"))], db=str(db), copy_to=str(scratch))
        return seen[0][1]

    def test_same_filesystem_not_copied(self, tmp_path, scratch, db, monkeypatch, caplog):
        monkeypatch.setattr(fcs_gx_purge, "_same_filesystem", lambda path, folder: True)
        assert self._db_used(tmp_path, db, scratch) == str(db)
        assert "same filesystem" in caplog.text and not scratch.exists()

    def test_network_copy_to_warned_but_copied(self, tmp_path, scratch, db, caplog):
        # the fixture reports every path as lustre
        assert Path(self._db_used(tmp_path, db, scratch)).parent.parent == scratch
        assert "network filesystem (lustre)" in caplog.text

    def test_same_filesystem_real_check(self, tmp_path):
        from aaftf.fcs_gx_purge import _same_filesystem

        f = tmp_path / "f"
        f.write_text("")
        assert _same_filesystem(str(f), str(tmp_path / "not" / "yet"))
