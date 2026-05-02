"""Tests for the AAFTF CLI framework.

Covers ALIAS_MAP routing, subcommand registration, argparse defaults,
required-argument enforcement, and the run_subtool() dispatcher.

No external bioinformatics tools are invoked.
"""

import sys
from argparse import Namespace
from unittest.mock import MagicMock, patch

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
            assert ALIAS_MAP.get(alias) == canonical, f"Alias '{alias}' should map to '{canonical}', " f"got '{ALIAS_MAP.get(alias)}'"

    def test_no_alias_maps_to_itself(self):
        for alias, canonical in ALIAS_MAP.items():
            assert alias != canonical, f"Alias '{alias}' maps to itself — should be a canonical name"


# ---------------------------------------------------------------------------
# run_subtool dispatcher
# ---------------------------------------------------------------------------


class TestRunSubtool:
    """Verify run_subtool() resolves aliases and handles unknown commands.

    We test ALIAS_MAP-level behaviour only.  The import-then-call dispatch
    cannot be reliably intercepted via sys.modules patching once submodules
    have been imported (Python binds 'import A.B as x' to the package
    attribute, not to sys.modules["A.B"]).  Full dispatch is exercised
    indirectly by the CLI parser tests that mock run_subtool at main().
    """

    def test_unknown_command_calls_parser_parse_args(self):
        """An unrecognised command should call parser.parse_args("")."""
        mock_parser = MagicMock()
        run_subtool(mock_parser, Namespace(command="definitely_not_real_xyzzy"))
        mock_parser.parse_args.assert_called_with("")

    def test_alias_map_resolves_stats_to_assess(self):
        assert ALIAS_MAP.get("stats") == "assess"

    def test_alias_map_resolves_coverage_to_depth(self):
        assert ALIAS_MAP.get("coverage") == "depth"

    def test_alias_map_resolves_cov_to_depth(self):
        assert ALIAS_MAP.get("cov") == "depth"

    def test_alias_map_resolves_dedup_to_rmdup(self):
        assert ALIAS_MAP.get("dedup") == "rmdup"


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
        args = _parse_with_main(["AAFTF", "depth", "-i", "genome.fa", "-l", "left.fq"])
        assert args.input == "genome.fa"

    def test_parses_left_reads(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "left.fq"])
        assert args.left == "left.fq"

    def test_parses_right_reads(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-r", "r.fq"])
        assert args.right == "r.fq"

    def test_parses_longreads(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq"])
        assert args.longreads == "lr.fq"

    def test_default_out(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.out == "coverage_stats.txt"

    def test_custom_out(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-o", "my_report.txt"])
        assert args.out == "my_report.txt"

    def test_default_cpus(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-c", "8"])
        assert args.cpus == 8

    def test_default_aligner(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.aligner == "minimap2"

    def test_bwa_aligner(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--aligner", "bwa"])
        assert args.aligner == "bwa"

    def test_default_longread_preset(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq"])
        assert args.longread_preset == "map-ont"

    def test_longread_preset_map_pb(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq", "--longread_type", "map-pb"])
        assert args.longread_preset == "map-pb"

    def test_debug_flag_false_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.debug is False

    def test_debug_flag_set(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-v"])
        assert args.debug is True

    def test_pipe_flag_false_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.pipe is False

    def test_missing_input_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "depth", "-l", "l.fq"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_coverage_alias_parses(self):
        args = _parse_with_main(["AAFTF", "coverage", "-i", "g.fa", "-l", "l.fq"])
        # alias 'coverage' routes to depth; same Namespace structure
        assert args.input == "g.fa"

    def test_default_plot_format_is_pdf(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.plot_format == "pdf"

    def test_plot_format_svg(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--plot-format", "svg"])
        assert args.plot_format == "svg"

    def test_plot_format_png(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--plot-format", "png"])
        assert args.plot_format == "png"

    def test_no_plot_false_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.no_plot is False

    def test_no_plot_flag(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--no-plot"])
        assert args.no_plot is True

    def test_default_quantize(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.quantize == "0:1:4:100:200:"

    def test_custom_quantize(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--quantize", "0:5:50:"])
        assert args.quantize == "0:5:50:"

    def test_quantize_labels_none_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.quantize_labels is None

    def test_custom_quantize_labels(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--quantize-labels", "NONE,LOW,HIGH"])
        assert args.quantize_labels == "NONE,LOW,HIGH"


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
        args = _parse_with_main(["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "-v"])
        assert args.debug is True

    def test_has_pipe_flag(self):
        args = _parse_with_main(["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "--pipe"])
        assert args.pipe is True

    def test_default_name_prefix(self):
        args = _parse_with_main(["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa"])
        assert args.name == "scaffold"


# ---------------------------------------------------------------------------
# fix_tbl parser — common flags present (added in audit 2026-05-02)
# ---------------------------------------------------------------------------


class TestFixTblParser:
    def test_has_debug_flag(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out), "-v"])
        assert args.debug is True

    def test_has_pipe_flag(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out), "--pipe"])
        assert args.pipe is True

    def test_debug_false_by_default(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out)])
        assert args.debug is False

    def test_pipe_false_by_default(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out)])
        assert args.pipe is False


# ---------------------------------------------------------------------------
# assess parser — telomere_window argument (added in audit 2026-05-02)
# ---------------------------------------------------------------------------


class TestAssessWindowParser:
    def test_default_telomere_window(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa"])
        assert args.telomere_window == 200

    def test_custom_telomere_window(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa", "--telomere_window", "500"])
        assert args.telomere_window == 500
