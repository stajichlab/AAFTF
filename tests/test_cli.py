"""Tests for the AAFTF CLI framework.

Covers ALIAS_MAP routing, subcommand registration, argparse defaults,
required-argument enforcement, and the run_subtool() dispatcher.

No external bioinformatics tools are invoked.
"""

import sys
from argparse import Namespace
from unittest.mock import patch, MagicMock

import pytest

from AAFTF.AAFTF_main import ALIAS_MAP, main, run_subtool

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# ALIAS_MAP
# ---------------------------------------------------------------------------

EXPECTED_ALIASES = {
    # trim
    "trim_reads": "trim",
    "read_trim": "trim",
    # filter
    "filter_reads": "filter",
    "read_filter": "filter",
    # assemble
    "asm": "assemble",
    "spades": "assemble",
    # vecscreen
    "vectorscreen": "vecscreen",
    "vector_blast": "vecscreen",
    # fcs_screen
    "ncbi_fcs": "fcs_screen",
    "ncbi_fcs-screen": "fcs_screen",
    # fcs_gx_purge
    "ncbi_fcs-gx": "fcs_gx_purge",
    "ncbi_fcs_gx": "fcs_gx_purge",
    "gx": "fcs_gx_purge",
    # sourpurge
    "purge": "sourpurge",
    # rmdup
    "dedup": "rmdup",
    # polish
    "pilon": "polish",
    "polca": "polish",
    # assess
    "stats": "assess",
    # fix_tbl
    "fix": "fix_tbl",
    # mito
    "mito_asm": "mito",
    "mitochondria": "mito",
    # depth
    "coverage": "depth",
    "cov": "depth",
}


class TestAliasMap:
    def test_depth_alias_coverage(self):
        assert ALIAS_MAP.get("coverage") == "depth"

    def test_depth_alias_cov(self):
        assert ALIAS_MAP.get("cov") == "depth"

    def test_all_expected_aliases_present(self):
        for alias, canonical in EXPECTED_ALIASES.items():
            assert ALIAS_MAP.get(alias) == canonical, (
                f"Alias '{alias}' should map to '{canonical}', "
                f"got '{ALIAS_MAP.get(alias)}'"
            )

    def test_no_alias_maps_to_itself(self):
        for alias, canonical in ALIAS_MAP.items():
            assert alias != canonical, (
                f"Alias '{alias}' maps to itself — should be a canonical name"
            )


# ---------------------------------------------------------------------------
# run_subtool dispatcher
# ---------------------------------------------------------------------------

class TestRunSubtool:
    """run_subtool() should import and call the right submodule.run()."""

    def _run_with_mock(self, command, module_path):
        mock_mod = MagicMock()
        args = Namespace(command=command)
        with patch.dict("sys.modules", {module_path: mock_mod}):
            run_subtool(None, args)
        return mock_mod

    def test_alias_routed_via_alias_map(self):
        """'stats' is an alias for 'assess'; run_subtool should load assess."""
        mock_mod = MagicMock()
        args = Namespace(command="stats")
        with patch("AAFTF.AAFTF_main.__import__", create=True):
            with patch.dict("sys.modules", {"AAFTF.assess": mock_mod}):
                run_subtool(None, args)
        mock_mod.run.assert_called_once_with(None, args)

    def test_depth_command_dispatched(self):
        mock_mod = MagicMock()
        args = Namespace(command="depth")
        with patch.dict("sys.modules", {"AAFTF.depth": mock_mod}):
            run_subtool(None, args)
        mock_mod.run.assert_called_once_with(None, args)

    def test_coverage_alias_dispatched_to_depth(self):
        mock_mod = MagicMock()
        args = Namespace(command="coverage")   # alias
        with patch.dict("sys.modules", {"AAFTF.depth": mock_mod}):
            run_subtool(None, args)
        mock_mod.run.assert_called_once_with(None, args)


# ---------------------------------------------------------------------------
# Subcommand help (argparse –– no tool execution)
# ---------------------------------------------------------------------------

def _parse_with_main(argv):
    """Run main() with sys.argv patched, intercept run_subtool before it fires."""
    captured = {}

    def _capture(parser, args):
        captured["args"] = args

    with patch.object(sys, "argv", argv):
        with patch("AAFTF.AAFTF_main.run_subtool", side_effect=_capture):
            main()
    return captured.get("args")


class TestSubcommandRegistration:
    def test_top_level_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0

    def test_depth_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "depth", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0

    def test_assess_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "assess", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0

    def test_sort_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "sort", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0


# ---------------------------------------------------------------------------
# depth parser — argument parsing and defaults
# ---------------------------------------------------------------------------

class TestDepthParser:
    def test_parses_input_flag(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "genome.fa", "-l", "left.fq"]
        )
        assert args.input == "genome.fa"

    def test_parses_left_reads(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "left.fq"]
        )
        assert args.left == "left.fq"

    def test_parses_right_reads(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-r", "r.fq"]
        )
        assert args.right == "r.fq"

    def test_parses_longreads(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq"]
        )
        assert args.longreads == "lr.fq"

    def test_default_out(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"]
        )
        assert args.out == "coverage_stats.txt"

    def test_custom_out(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-o", "my_report.txt"]
        )
        assert args.out == "my_report.txt"

    def test_default_cpus(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"]
        )
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-c", "8"]
        )
        assert args.cpus == 8

    def test_default_aligner(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"]
        )
        assert args.aligner == "minimap2"

    def test_bwa_aligner(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--aligner", "bwa"]
        )
        assert args.aligner == "bwa"

    def test_default_longread_preset(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq"]
        )
        assert args.longread_preset == "map-ont"

    def test_longread_preset_map_pb(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq",
             "--longread_type", "map-pb"]
        )
        assert args.longread_preset == "map-pb"

    def test_debug_flag_false_by_default(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"]
        )
        assert args.debug is False

    def test_debug_flag_set(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-v"]
        )
        assert args.debug is True

    def test_pipe_flag_false_by_default(self):
        args = _parse_with_main(
            ["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"]
        )
        assert args.pipe is False

    def test_missing_input_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "depth", "-l", "l.fq"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_coverage_alias_parses(self):
        args = _parse_with_main(
            ["AAFTF", "coverage", "-i", "g.fa", "-l", "l.fq"]
        )
        # alias 'coverage' routes to depth; same Namespace structure
        assert args.input == "g.fa"


# ---------------------------------------------------------------------------
# assess parser — common flags present
# ---------------------------------------------------------------------------

class TestAssessParser:
    def test_has_debug_flag(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa", "-v"])
        assert args.debug is True

    def test_has_pipe_flag(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa", "--pipe"])
        assert args.pipe is True

    def test_default_telomere_monomer(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa"])
        assert args.telomere_monomer == "TAA[C]+"

    def test_default_telomere_n_repeat(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa"])
        assert args.telomere_n_repeat == 2


# ---------------------------------------------------------------------------
# sort parser — common flags present
# ---------------------------------------------------------------------------

class TestSortParser:
    def test_has_debug_flag(self):
        args = _parse_with_main(
            ["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "-v"]
        )
        assert args.debug is True

    def test_has_pipe_flag(self):
        args = _parse_with_main(
            ["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "--pipe"]
        )
        assert args.pipe is True

    def test_default_name_prefix(self):
        args = _parse_with_main(
            ["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa"]
        )
        assert args.name == "scaffold"
