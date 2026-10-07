"""Unit tests for aaftf/dependency.py."""

import pytest

from aaftf.dependency import _print_tool_table

pytestmark = pytest.mark.unit

RESULTS = [("bwa", "/opt/bin/bwa", "filter, depth"), ("bowtie2", None, "filter")]


def test_full_table_shows_location_and_users(capsys):
    assert _print_tool_table("Checking...", RESULTS) == ["bowtie2"]
    out = capsys.readouterr().out
    assert "/opt/bin/bwa" in out and "used by: filter" in out


def test_quiet_table_shows_only_status_and_name(capsys):
    assert _print_tool_table("Checking...", RESULTS, quiet=True) == ["bowtie2"]
    assert capsys.readouterr().out.splitlines() == ["  [OK]      bwa", "  [MISSING] bowtie2"]


class TestRun:
    """Only tools and packages the default pipeline needs make `dependency` fail."""

    def _run(self, missing_tools=(), missing_imports=()):
        from unittest.mock import patch

        from aaftf import dependency

        real_import = dependency.importlib.import_module

        def _which(tool):
            return None if tool in missing_tools else f"/opt/bin/{tool}"

        def _import(name):
            if name in missing_imports:
                raise ImportError(name)
            return real_import(name)

        with patch("aaftf.dependency.shutil.which", side_effect=_which), patch("aaftf.dependency.importlib.import_module", side_effect=_import):
            dependency.run()

    def test_everything_installed(self, caplog, capsys):
        import logging

        with caplog.at_level(logging.INFO, logger="aaftf.dependency"):
            self._run()
        assert "Everything the default pipeline needs is installed" in caplog.text
        assert "optional tool(s)/package(s) missing" not in caplog.text
        assert "[MISSING]" not in capsys.readouterr().out

    def test_missing_optional_tool_and_package_are_not_errors(self, caplog, capsys):
        import logging

        with caplog.at_level(logging.INFO, logger="aaftf.dependency"):
            self._run(missing_tools={"bowtie2", "megahit", "mosdepth"}, missing_imports={"matplotlib"})
        assert "Everything the default pipeline needs is installed" in caplog.text
        (note,) = (r.getMessage() for r in caplog.records if "optional tool(s)/package(s) missing" in r.getMessage())
        assert note.startswith("NOTE: 4 ")
        for name in ("bowtie2", "megahit", "mosdepth", "matplotlib"):
            assert name in note
        out = capsys.readouterr().out
        assert "[MISSING] bowtie2" in out and "[MISSING] matplotlib" in out

    @pytest.mark.parametrize("tool", ["spades.py", "bbduk.sh", "bwa"])
    def test_missing_pipeline_tool_is_an_error(self, tool):
        with pytest.raises(FileNotFoundError, match=tool):
            self._run(missing_tools={tool})

    def test_missing_required_package_is_an_error(self):
        with pytest.raises(FileNotFoundError, match="psutil"):
            self._run(missing_imports={"psutil"})


def test_required_and_optional_tools_disjoint():
    from aaftf.dependency import OPTIONAL_TOOLS, REQUIRED_TOOLS

    assert not set(REQUIRED_TOOLS) & set(OPTIONAL_TOOLS)


def test_tools_checked_in_written_order():
    from unittest.mock import patch

    from aaftf.dependency import OPTIONAL_TOOLS, _check_tools

    with patch("aaftf.dependency.shutil.which", return_value=None):
        assert [tool for tool, _, _ in _check_tools(OPTIONAL_TOOLS)] == list(OPTIONAL_TOOLS)


@pytest.mark.parametrize("avx2_binary, suggested", [(False, True), (True, False)])
def test_bowtie2_avx2_suggestion(tmp_path, caplog, avx2_binary, suggested):
    """The install-bowtie2 hint appears only for the conda bowtie2 (no bowtie2-align-s-v256 next to it)."""
    import logging
    from unittest.mock import patch

    from aaftf.dependency import _suggest_bowtie2_avx2

    (tmp_path / "bowtie2").write_text("")
    if avx2_binary:
        (tmp_path / "bowtie2-align-s-v256").write_text("")
    with patch("aaftf.dependency.shutil.which", return_value=str(tmp_path / "bowtie2")), caplog.at_level(logging.INFO, logger="aaftf"):
        _suggest_bowtie2_avx2()
    assert ("install-bowtie2" in caplog.text) is suggested
