"""Unit tests for aaftf/fix_tbl.py.

All functions are pure Python — no external tools required.
"""

import io
import logging

import pytest

from aaftf.fix_tbl import TblFeature, fix_tbl, parse_adjustments, parse_tbl

pytestmark = pytest.mark.unit


def _fh(text):
    """Wrap *text* in a StringIO for use as an input file handle."""
    return io.StringIO(text)


def _report(*rows):
    """Build an FCS action report with the header line followed by tab-joined *rows*."""
    lines = ['##["FCS-adaptor", 1, 1]', "#accession\tlength\taction\trange\tname"]
    return _fh("\n".join(lines + ["\t".join(map(str, row)) for row in rows]) + "\n")


def _fix(tbl, *rows):
    """Run fix_tbl on *tbl* with an FCS report of *rows*; return {seqid: [TblFeature]} parsed from the output."""
    out = io.StringIO()
    fix_tbl(_fh(tbl), _report(*rows), out)
    return parse_tbl(_fh(out.getvalue()))


# ---------------------------------------------------------------------------
# parse_tbl
# ---------------------------------------------------------------------------

TBL_BASIC = """\
>Feature scaffold_1
1\t1000\tgene
\t\t\tlocus_tag\tGENE_001
1\t300\tCDS
400\t900
\t\t\tproduct\thypothetical protein
>Feature scaffold_2
1500\t500\tgene
\t\t\tlocus_tag\tGENE_002
"""


class TestParseTbl:
    def test_sequences_in_file_order(self):
        assert list(parse_tbl(_fh(TBL_BASIC))) == ["scaffold_1", "scaffold_2"]

    def test_features_group_intervals_and_qualifiers(self):
        gene, cds = parse_tbl(_fh(TBL_BASIC))["scaffold_1"]
        assert gene == TblFeature("gene", [["1", "1000"]], ["\t\t\tlocus_tag\tGENE_001"])
        assert cds.intervals == [["1", "300"], ["400", "900"]]
        assert cds.qualifiers == ["\t\t\tproduct\thypothetical protein"]

    def test_minus_strand_kept_as_written(self):
        assert parse_tbl(_fh(TBL_BASIC))["scaffold_2"][0].intervals == [["1500", "500"]]

    def test_blank_and_comment_lines_skipped(self):
        features = parse_tbl(_fh("# comment\n>Feature s1\n\n1\t100\tgene\n\n>Feature s2\n"))
        assert len(features["s1"]) == 1 and features["s2"] == []

    def test_partial_markers_kept(self):
        feat = parse_tbl(_fh(">Feature s1\n<1\t>1000\tgene\n"))["s1"][0]
        assert feat.intervals == [["<1", ">1000"]]

    def test_feature_before_header_raises(self):
        with pytest.raises(ValueError, match="before the first '>Feature'"):
            parse_tbl(_fh("1\t100\tgene\n>Feature s1\n"))

    def test_round_trip(self):
        features = parse_tbl(_fh(TBL_BASIC))
        text = "".join(f">Feature {seqid}\n" + "".join(f.to_tbl() for f in feats) for seqid, feats in features.items())
        assert text == TBL_BASIC


# ---------------------------------------------------------------------------
# parse_adjustments
# ---------------------------------------------------------------------------


class TestParseAdjustments:
    def test_first_row_after_header_is_read(self):
        adj = parse_adjustments(_report(("scaffold_1", 1000, "ACTION_TRIM", "1..50", "adaptor")))
        assert adj == {"scaffold_1": [(1000, "ACTION_TRIM", 1, 50)]}

    def test_multiple_ranges(self):
        adj = parse_adjustments(_report(("scaffold_1", 1000, "ACTION_TRIM", "1..30,971..1000", "adaptor")))
        assert adj["scaffold_1"] == [(1000, "ACTION_TRIM", 1, 30), (1000, "ACTION_TRIM", 971, 1000)]

    def test_exclude_without_range_covers_whole_sequence(self):
        adj = parse_adjustments(_report(("scaffold_1", 1000, "ACTION_EXCLUDE", "", "adaptor")))
        assert adj["scaffold_1"] == [(1000, "ACTION_EXCLUDE", 1, 1000)]

    def test_lines_before_header_ignored(self):
        assert parse_adjustments(_fh("scaffold_1\t1000\tACTION_TRIM\t1..50\n")) == {}


# ---------------------------------------------------------------------------
# fix_tbl
# ---------------------------------------------------------------------------

TBL_FOR_FIX = """\
>Feature scaffold_1
10\t40\tgene
\t\t\tlocus_tag\tGENE_IN_TRIM
30\t100\tgene
\t\t\tlocus_tag\tGENE_CUT_LEFT
200\t300\tgene
\t\t\tlocus_tag\tGENE_B
300\t200\tgene
\t\t\tlocus_tag\tGENE_MINUS
750\t900\tgene
\t\t\tlocus_tag\tGENE_CUT_RIGHT
850\t950\tgene
\t\t\tlocus_tag\tGENE_PAST_END
>Feature scaffold_2
1\t100\tgene
"""


def _by_tag(features, seqid="scaffold_1"):
    """Map locus_tag → intervals for the features of *seqid*."""
    return {f.qualifiers[0].split("\t")[-1]: f.intervals for f in features[seqid]}


class TestFixTblLeftTrim:
    ROW = ("scaffold_1", 1000, "ACTION_TRIM", "1..50", "adaptor")

    def test_features_shifted(self):
        assert _by_tag(_fix(TBL_FOR_FIX, self.ROW))["GENE_B"] == [["150", "250"]]

    def test_minus_strand_shifted(self):
        assert _by_tag(_fix(TBL_FOR_FIX, self.ROW))["GENE_MINUS"] == [["250", "150"]]

    def test_feature_cut_by_trim_marked_5prime_partial(self):
        assert _by_tag(_fix(TBL_FOR_FIX, self.ROW))["GENE_CUT_LEFT"] == [["<1", "50"]]

    def test_feature_inside_trim_dropped(self, caplog):
        with caplog.at_level(logging.WARNING, logger="aaftf.fix_tbl"):
            tags = _by_tag(_fix(TBL_FOR_FIX, self.ROW))
        assert "GENE_IN_TRIM" not in tags
        assert "inside a trimmed region" in caplog.text

    def test_multiple_left_ranges_use_the_furthest(self):
        rows = [("scaffold_1", 1000, "ACTION_TRIM", "1..20", "a"), self.ROW]
        assert _by_tag(_fix(TBL_FOR_FIX, *rows))["GENE_B"] == [["150", "250"]]


class TestFixTblRightTrim:
    ROW = ("scaffold_1", 1000, "ACTION_TRIM", "801..1000", "adaptor")

    def test_features_before_trim_unchanged(self):
        assert _by_tag(_fix(TBL_FOR_FIX, self.ROW))["GENE_B"] == [["200", "300"]]

    def test_feature_cut_by_trim_clipped_and_marked_3prime_partial(self):
        assert _by_tag(_fix(TBL_FOR_FIX, self.ROW))["GENE_CUT_RIGHT"] == [["750", ">800"]]

    def test_feature_past_new_end_dropped(self):
        assert "GENE_PAST_END" not in _by_tag(_fix(TBL_FOR_FIX, self.ROW))

    def test_minus_strand_cut_marks_5prime_end(self):
        tbl = ">Feature s1\n900\t700\tgene\n\t\t\tlocus_tag\tM\n"
        assert _by_tag(_fix(tbl, ("s1", 1000, "ACTION_TRIM", "801..1000", "a")), "s1")["M"] == [["<800", "700"]]


class TestFixTblMultiInterval:
    TBL = ">Feature s1\n20\t100\tCDS\n200\t300\n400\t500\n\t\t\tlocus_tag\tX\n"

    def test_dropped_first_interval_marks_5prime_partial(self):
        assert _by_tag(_fix(self.TBL, ("s1", 1000, "ACTION_TRIM", "1..150", "a")), "s1")["X"] == [["<50", "150"], ["250", "350"]]

    def test_dropped_last_interval_marks_3prime_partial(self):
        assert _by_tag(_fix(self.TBL, ("s1", 1000, "ACTION_TRIM", "351..1000", "a")), "s1")["X"] == [["20", "100"], ["200", ">300"]]


class TestFixTblOther:
    def test_existing_partial_markers_kept(self):
        tbl = ">Feature s1\n<100\t>200\tgene\n\t\t\tlocus_tag\tP\n"
        assert _by_tag(_fix(tbl, ("s1", 1000, "ACTION_TRIM", "1..50", "a")), "s1")["P"] == [["<50", ">150"]]

    def test_untrimmed_sequence_unchanged(self):
        row = ("scaffold_1", 1000, "ACTION_TRIM", "1..50", "adaptor")
        assert _fix(TBL_FOR_FIX, row)["scaffold_2"][0].intervals == [["1", "100"]]

    def test_excluded_sequence_dropped(self, caplog):
        with caplog.at_level(logging.WARNING, logger="aaftf.fix_tbl"):
            features = _fix(TBL_FOR_FIX, ("scaffold_2", 1000, "ACTION_EXCLUDE", "", "adaptor"))
        assert "scaffold_2" not in features
        assert "excluded" in caplog.text

    def test_internal_trim_left_alone_with_warning(self, caplog):
        with caplog.at_level(logging.WARNING, logger="aaftf.fix_tbl"):
            features = _fix(TBL_FOR_FIX, ("scaffold_1", 1000, "ACTION_TRIM", "400..500", "adaptor"))
        assert _by_tag(features)["GENE_B"] == [["200", "300"]]
        assert "inside the contig" in caplog.text

    def test_qualifiers_preserved(self):
        out = io.StringIO()
        fix_tbl(_fh(TBL_FOR_FIX), _report(), out)
        assert out.getvalue() == TBL_FOR_FIX


# ---------------------------------------------------------------------------
# run (file paths in, file path out)
# ---------------------------------------------------------------------------


class TestRun:
    def _files(self, tmp_path, tbl_text):
        tbl, report, out = tmp_path / "in.tbl", tmp_path / "fcs_report.txt", tmp_path / "out.tbl"
        tbl.write_text(tbl_text)
        report.write_text("#accession\tlength\taction\trange\tname\nscaffold_1\t1000\tACTION_TRIM\t1..50\tadaptor\n")
        return tbl, report, out

    def test_writes_fixed_table(self, tmp_path):
        from aaftf.fix_tbl import run

        tbl, report, out = self._files(tmp_path, ">Feature scaffold_1\n200\t300\tgene\n\t\t\tlocus_tag\tA\n")
        run(table=str(tbl), report=str(report), output=str(out))
        assert out.read_text() == ">Feature scaffold_1\n150\t250\tgene\n\t\t\tlocus_tag\tA\n"
        assert not (tmp_path / "out.tbl.tmp").exists()

    def test_malformed_input_leaves_existing_output_untouched(self, tmp_path):
        from aaftf.fix_tbl import run

        tbl, report, out = self._files(tmp_path, "1\t100\tgene\n")  # feature line before any header
        out.write_text("previous result\n")
        with pytest.raises(ValueError):
            run(table=str(tbl), report=str(report), output=str(out))
        assert out.read_text() == "previous result\n"
        assert not (tmp_path / "out.tbl.tmp").exists()
