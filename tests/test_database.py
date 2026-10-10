"""Unit tests for aaftf/database.py (no network access)."""

import os
import stat
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.database import _resolve, run
from aaftf.main import main
from aaftf.resources import DATABASES

pytestmark = pytest.mark.unit


def _row(out, name):
    """Return the listing row whose first column is ``name``."""
    return next(line for line in out.splitlines() if line.split() and line.split()[0] == name)


def _fake_download(calls):
    """Stand-in for download_file that records calls and writes a stub file."""

    def _download(url, dest, force=False):
        calls.append((url, dest))
        if force or not Path(dest).exists():
            Path(dest).parent.mkdir(parents=True, exist_ok=True)
            Path(dest).write_text("stub")
        return dest

    return _download


class TestListDb:
    def _list(self, monkeypatch, capsys, folders):
        monkeypatch.setenv("AAFTF_DB", os.pathsep.join(str(f) for f in folders))
        with patch("aaftf.database._remote_size", return_value=1024):
            run()
        return capsys.readouterr().out

    def test_no_arguments_lists_every_database(self, monkeypatch, capsys, tmp_path):
        out = self._list(monkeypatch, capsys, [tmp_path])
        for abbr, entry in DATABASES.items():
            row = _row(out, abbr)
            assert entry["filename"] in row and entry["used_by"] in row

    def test_shows_folder_each_file_is_stored_in(self, monkeypatch, capsys, tmp_path):
        shared, mine = tmp_path / "shared", tmp_path / "mine"
        shared.mkdir()
        mine.mkdir()
        (shared / DATABASES["univec"]["filename"]).write_text("x")
        (mine / DATABASES["proks"]["filename"]).write_text("x")

        out = self._list(monkeypatch, capsys, [shared, mine])

        assert _row(out, "univec").endswith(str(shared.resolve()))
        assert _row(out, "proks").endswith(str(mine.resolve()))
        assert _row(out, "euks").endswith("not downloaded")
        assert "1.0 KB" in _row(out, "euks")  # remote size estimate

    def test_first_folder_wins_for_duplicates(self, monkeypatch, capsys, tmp_path):
        shared, mine = tmp_path / "shared", tmp_path / "mine"
        for folder in (shared, mine):
            folder.mkdir()
            (folder / DATABASES["univec"]["filename"]).write_text("x")

        assert _row(self._list(monkeypatch, capsys, [shared, mine]), "univec").endswith(str(shared.resolve()))

    def test_other_files_show_their_folder(self, monkeypatch, capsys, tmp_path):
        (tmp_path / "notes.txt").write_text("hello")
        row = next(line for line in self._list(monkeypatch, capsys, [tmp_path]).splitlines() if "notes.txt" in line)
        assert row.endswith(f"{tmp_path.resolve()} (other file)")

    def test_list_does_not_create_folders(self, monkeypatch, capsys, tmp_path):
        missing = tmp_path / "missing"
        self._list(monkeypatch, capsys, [missing])
        assert not missing.exists()


class TestResolve:
    def test_abbreviation_and_file_name(self):
        assert _resolve(["univec", "contam_in_prok.fa"]) == ["univec", "proks"]

    def test_case_insensitive(self):
        assert _resolve(["UNIVEC", "Contam_In_Prok.FA"]) == ["univec", "proks"]

    def test_all_and_duplicates(self):
        selected = _resolve(["univec", "all", "UniVec"])
        assert selected[0] == "univec" and sorted(selected) == sorted(DATABASES) and len(selected) == len(set(selected))

    def test_sourmash_database_brings_its_lineages(self):
        assert _resolve(["gtdb-reps-lite"]) == ["gtdb-reps-lite", "gtdb-reps-lineages"]
        assert _resolve(["gtdb-reps", "gtdb-reps-lite"]) == ["gtdb-reps", "gtdb-reps-lineages", "gtdb-reps-lite"]
        assert _resolve(["human-hg38"]) == ["human-hg38"]  # taxpurge writes its lineage itself

    def test_taxpurge_keyword_is_its_default_databases(self):
        assert _resolve(["taxpurge"]) == ["gbk-fungi-lite", "gbk-euks-lineages", "gtdb-reps-lite", "gtdb-reps-lineages", "human-hg38"]

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError) as exc:
            _resolve(["univec", "not_a_db"])
        assert "not_a_db" in str(exc.value)


class TestDownload:
    def test_downloads_named_databases_to_write_folder(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        calls = []
        with patch("aaftf.database.download_file", side_effect=_fake_download(calls)):
            run(databases=["univec", "contam_in_prok.fa"])
        assert calls == [
            (DATABASES["univec"]["url"], str(tmp_path.resolve() / "UniVec")),
            (DATABASES["proks"]["url"], str(tmp_path.resolve() / "contam_in_prok.fa")),
        ]

    def test_existing_copy_in_another_folder_is_reused(self, monkeypatch, tmp_path):
        shared, mine = tmp_path / "shared", tmp_path / "mine"
        shared.mkdir()
        (shared / "UniVec").write_text("shared")
        monkeypatch.setenv("AAFTF_DB", f"{shared}{os.pathsep}{mine}")
        calls = []
        with patch("aaftf.database.download_file", side_effect=_fake_download(calls)):
            run(databases=["univec"])
        assert calls[0][1] == str(shared.resolve() / "UniVec")

    def test_fcs_script_made_executable(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        with patch("aaftf.database.download_file", side_effect=_fake_download([])):
            run(databases=["fcs_script"])
        assert (tmp_path / "run_fcsadaptor.sh").stat().st_mode & stat.S_IXUSR

    def test_failed_download_raises_after_trying_the_rest(self, monkeypatch, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        calls = []
        ok = _fake_download(calls)

        def _flaky(url, dest, force=False):
            if url == DATABASES["univec"]["url"]:
                raise OSError("network down")
            return ok(url, dest, force)

        with patch("aaftf.database.download_file", side_effect=_flaky):
            with pytest.raises(RuntimeError):
                run(databases=["univec", "proks"])
        assert [url for url, _ in calls] == [DATABASES["proks"]["url"]]


class TestDownloadCli:
    def _parse(self, argv):
        captured = {}
        with patch.object(sys, "argv", argv):
            with patch("aaftf.database.run", side_effect=lambda **kw: captured.update(args=Namespace(**kw))):
                main()
        return captured["args"]

    def test_no_arguments_means_list(self):
        assert self._parse(["AAFTF", "database"]).databases == []

    def test_names_and_force(self):
        args = self._parse(["AAFTF", "database", "univec", "contam_in_prok.fa", "--force"])
        assert args.databases == ["univec", "contam_in_prok.fa"] and args.force is True

    def test_old_flags_removed(self):
        with patch.object(sys, "argv", ["AAFTF", "database", "--skip-core"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0


class TestRequiredOptionalGroups:
    """Databases are split, like `AAFTF dependency`, into what the default pipeline needs and the rest."""

    def test_every_entry_is_marked(self):
        assert all(isinstance(entry["required"], bool) for entry in DATABASES.values())

    def test_required_set_is_what_the_default_pipeline_reads(self):
        # filter reads phix + univec; vecscreen univec, euks, proks, mitodb (taxpurge and its sourmash databases run only with -p)
        assert {abbr for abbr, e in DATABASES.items() if e["required"]} == {"phix", "univec", "euks", "proks", "mitodb"}

    def test_listing_shows_required_then_optional(self, monkeypatch, capsys, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        with patch("aaftf.database._remote_size", return_value=1024):
            run()
        out = capsys.readouterr().out
        required_at, optional_at = out.index("Needed by the default pipeline:"), out.index("Only for non-default options or optional steps:")
        assert required_at < out.index("\n  phix ") < out.index("\n  mitodb ") < optional_at < out.index("\n  gbk-fungi ") < out.index("\n  human-hg38 ") < out.index("\n  fcs_image ")
        assert "Other files" not in out  # heading only when the folder holds other files

    def test_other_files_get_their_own_heading(self, monkeypatch, capsys, tmp_path):
        (tmp_path / "notes.txt").write_text("x")
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        with patch("aaftf.database._remote_size", return_value=1024):
            run()
        out = capsys.readouterr().out
        assert out.index("Other files in the database folders:") < out.index("notes.txt")

    def _listing(self, monkeypatch, capsys, db):
        monkeypatch.setenv("AAFTF_DB", str(db))
        with patch("aaftf.database._remote_size", return_value=1024):
            run()
        return capsys.readouterr().out.splitlines()

    def _index(self, folder):
        folder.mkdir(parents=True)
        (folder / "CURRENT").write_text("MANIFEST-000005\n")

    def test_built_index_listed_under_its_sketch_file(self, monkeypatch, capsys, tmp_path):
        """The index row follows its sketch file's row; a sketch file replaced by it is "fulfilled by rocksdb"."""
        self._index(tmp_path / "hg38-entire.rocksdb")
        lines = self._listing(monkeypatch, capsys, tmp_path)
        row = next(i for i, line in enumerate(lines) if line.startswith("  human-hg38 "))
        assert "hg38-entire.sig.zip" in lines[row] and "fulfilled by rocksdb" in lines[row]
        assert "1.0 KB" in lines[row]  # the sketch file's download size (mocked), not "-"
        assert "hg38-entire.rocksdb" in lines[row + 1] and str(tmp_path.resolve()) in lines[row + 1]
        fungi = next(i for i, line in enumerate(lines) if line.startswith("  gbk-fungi "))
        assert "not downloaded" in lines[fungi] and "genbank-20250408-fungi-k31.rocksdb" in lines[fungi + 1] and "not built" in lines[fungi + 1]
        assert not any("Other files" in line for line in lines)

    def test_published_index_is_the_database_row(self, monkeypatch, capsys, tmp_path):
        self._index(tmp_path / "gtdb-rs226-reps.k31.rocksdb")
        lines = self._listing(monkeypatch, capsys, tmp_path)
        row = next(line for line in lines if line.startswith("  gtdb-reps "))
        assert "gtdb-rs226-reps.k31.rocksdb" in row and str(tmp_path.resolve()) in row
        assert not any(".sig.zip" in line and "gtdb-rs226-reps.k31" in line for line in lines)  # no sketch file entry

    def test_other_files_lists_only_the_top_level(self, monkeypatch, capsys, tmp_path):
        """A subfolder is one entry (name/); the files inside it are not listed."""
        (tmp_path / "old_run").mkdir()
        (tmp_path / "old_run" / "inner.sst").write_text("x")
        (tmp_path / "top.txt").write_text("x")
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        with patch("aaftf.database._remote_size", return_value=1024):
            run()
        out = capsys.readouterr().out
        assert "old_run/" in out and "top.txt" in out
        assert "inner.sst" not in out

    def test_required_and_optional_keywords(self):
        assert _resolve(["required"]) == ["phix", "univec", "euks", "proks", "mitodb"]
        assert _resolve(["optional"]) == [
            "gbk-fungi",
            "gbk-fungi-lineages",
            "gbk-fungi-lite",
            "gbk-euks-lineages",
            "gtdb-reps",
            "gtdb-reps-lineages",
            "gtdb-reps-lite",
            "human-hg38",
            "fcs_script",
            "fcs_image",
        ]
        assert _resolve(["REQUIRED", "fcs_script"]) == ["phix", "univec", "euks", "proks", "mitodb", "fcs_script"]


class TestInstallIndex:
    """Sourmash databases are installed as RocksDB indexes: unpacked from a published archive, or built."""

    def _archive(self, tmp_path, folder):
        import tarfile

        src = tmp_path / "src" / folder
        src.mkdir(parents=True)
        (src / "CURRENT").write_text("MANIFEST-000005\n")
        (src / "000001.sst").write_text("x")
        archive = tmp_path / f"{folder}.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(src, arcname=folder)
        return archive

    def _fake_download(self, archive=None):
        def download(url, dest, force=False):
            Path(dest).parent.mkdir(parents=True, exist_ok=True)
            Path(dest).write_bytes(archive.read_bytes() if archive and url.endswith(".tar.gz") else b"x")
            return dest

        return download

    def _fake_index(self, cmds):
        def index(cmd, debug=False, **kw):
            cmds.append(cmd)
            out = Path(cmd[cmd.index("-o") + 1])
            out.mkdir()
            (out / "CURRENT").write_text("MANIFEST-000005\n")

        return index

    def test_published_archive_unpacked_into_place(self, monkeypatch, tmp_path):
        db = tmp_path / "db"
        monkeypatch.setenv("AAFTF_DB", str(db))
        archive = self._archive(tmp_path, "gtdb-rs226-reps.k31-sc10k.sig.rocksdb")
        with patch("aaftf.database.download_file", side_effect=self._fake_download(archive)):
            run(databases=["gtdb-reps-lite"])
        index = db / "gtdb-rs226-reps.k31-sc10k.sig.rocksdb"
        assert (index / "CURRENT").is_file() and (index / "000001.sst").is_file()
        assert sorted(p.name for p in db.iterdir()) == ["gtdb-rs226-reps.k31-sc10k.sig.rocksdb", "gtdb-rs226-reps.lineages.csv"]  # archive and temp folder removed

    def test_archive_without_the_index_folder_fails_and_leaves_nothing(self, monkeypatch, tmp_path):
        db = tmp_path / "db"
        monkeypatch.setenv("AAFTF_DB", str(db))
        archive = self._archive(tmp_path, "something-else.rocksdb")
        with patch("aaftf.database.download_file", side_effect=self._fake_download(archive)):
            with pytest.raises(RuntimeError, match="gtdb-reps-lite"):
                run(databases=["gtdb-reps-lite"])
        assert not any(p.name.startswith("gtdb-rs226-reps.k31-sc10k.sig.rocksdb") for p in db.iterdir())

    def test_built_from_its_sketch_file_which_is_then_removed(self, monkeypatch, tmp_path):
        db = tmp_path / "db"
        monkeypatch.setenv("AAFTF_DB", str(db))
        cmds = []
        with patch("aaftf.database.download_file", side_effect=self._fake_download()), patch("aaftf.database.run_cmd", side_effect=self._fake_index(cmds)), patch("aaftf.database.require_tools"):
            run(databases=["human-hg38"])
        assert cmds[0][:4] == ["sourmash", "scripts", "index", str(db.resolve() / "hg38-entire.sig.zip")]
        assert cmds[0][-4:] == ["-k", "31", "-s", "1000"]
        assert (db / "hg38-entire.rocksdb" / "CURRENT").is_file()
        assert not (db / "hg38-entire.sig.zip").exists()  # replaced by the index

    def test_no_build_keeps_the_sketch_file(self, monkeypatch, tmp_path):
        db = tmp_path / "db"
        monkeypatch.setenv("AAFTF_DB", str(db))
        with patch("aaftf.database.download_file", side_effect=self._fake_download()), patch("aaftf.database.run_cmd") as rc:
            run(databases=["human-hg38"], no_build=True)
        rc.assert_not_called()
        assert (db / "hg38-entire.sig.zip").is_file() and not (db / "hg38-entire.rocksdb").exists()

    def test_sketch_file_in_another_folder_is_not_removed(self, monkeypatch, tmp_path):
        """A sketch file in a second (e.g. shared) $AAFTF_DB folder is used to build, but kept."""
        mine, shared = tmp_path / "mine", tmp_path / "shared"
        mine.mkdir()
        shared.mkdir()
        (shared / "hg38-entire.sig.zip").write_text("x")
        monkeypatch.setenv("AAFTF_DB", f"{mine}{os.pathsep}{shared}")
        with patch("aaftf.database.run_cmd", side_effect=self._fake_index([])), patch("aaftf.database.require_tools"):
            run(databases=["human-hg38"])
        assert (mine / "hg38-entire.rocksdb" / "CURRENT").is_file() and (shared / "hg38-entire.sig.zip").is_file()

    def test_already_present_index_is_skipped(self, monkeypatch, tmp_path):
        db = tmp_path / "db"
        (db / "hg38-entire.rocksdb").mkdir(parents=True)
        (db / "hg38-entire.rocksdb" / "CURRENT").write_text("x")
        monkeypatch.setenv("AAFTF_DB", str(db))
        with patch("aaftf.database.run_cmd") as rc, patch("aaftf.database.download_file") as dl:
            run(databases=["human-hg38"])
        rc.assert_not_called()
        dl.assert_not_called()
