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
        assert _resolve(["univec", "all", "UniVec"]) == ["univec"] + [a for a in DATABASES if a != "univec"]

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
        # filter reads phix + univec; vecscreen univec, euks, proks, mitodb (taxpurge and its sm_gbk run only with -p)
        assert {abbr for abbr, e in DATABASES.items() if e["required"]} == {"phix", "univec", "euks", "proks", "mitodb"}

    def test_listing_shows_required_then_optional(self, monkeypatch, capsys, tmp_path):
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        with patch("aaftf.database._remote_size", return_value=1024):
            run()
        out = capsys.readouterr().out
        required_at, optional_at = out.index("Needed by the default pipeline:"), out.index("Only for non-default options or optional steps:")
        assert required_at < out.index("\n  phix ") < out.index("\n  mitodb ") < optional_at < out.index("\n  sm_gbk ") < out.index("\n  sm_gtdb ") < out.index("\n  fcs_image ")
        assert "Other files" not in out  # heading only when the folder holds other files

    def test_other_files_get_their_own_heading(self, monkeypatch, capsys, tmp_path):
        (tmp_path / "notes.txt").write_text("x")
        monkeypatch.setenv("AAFTF_DB", str(tmp_path))
        with patch("aaftf.database._remote_size", return_value=1024):
            run()
        out = capsys.readouterr().out
        assert out.index("Other files in the database folders:") < out.index("notes.txt")

    def test_required_and_optional_keywords(self):
        assert _resolve(["required"]) == ["phix", "univec", "euks", "proks", "mitodb"]
        assert _resolve(["optional"]) == ["sm_gbk", "sm_gtdbrep", "sm_gtdb", "fcs_script", "fcs_image"]
        assert _resolve(["REQUIRED", "fcs_script"]) == ["phix", "univec", "euks", "proks", "mitodb", "fcs_script"]
