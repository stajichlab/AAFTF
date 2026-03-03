"""Unit tests for AAFTF/fix_tbl.py.

All functions are pure Python — no external tools required.
"""

import io

import pytest

from AAFTF.fix_tbl import fix_tbl, parse_adjustments, parse_tbl

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper builders
# ---------------------------------------------------------------------------

def _tbl(text):
    """Wrap *text* in a StringIO for use as a tbl file handle."""
    return io.StringIO(text)


def _adj(text):
    """Wrap *text* in a StringIO for use as an adjustment file handle."""
    return io.StringIO(text)


def _out():
    """Return a writable StringIO to capture output."""
    return io.StringIO()


# ---------------------------------------------------------------------------
# parse_tbl
# ---------------------------------------------------------------------------

TBL_BASIC = """\
>Feature scaffold_1
1\t1000\tgene
\t\t\tlocus_tag\tGENE_001
1\t900\tCDS
\t\t\tproduct\thypothetical protein
>Feature scaffold_2
500\t1500\tgene
\t\t\tlocus_tag\tGENE_002
"""

TBL_WITH_BLANKS = """\
>Feature scaffold_1

1\t500\tgene
\t\t\tlocus_tag\tGENE_003

>Feature scaffold_2
100\t200\tgene
"""

TBL_WITH_COMMENT = """\
# This is a comment
>Feature scaffold_1
1\t100\tgene
"""

TBL_PARTIAL_COORDS = """\
>Feature scaffold_1
<1\t>1000\tgene
\t\t\tlocus_tag\tGENE_005
"""


class TestParseTbl:
    def test_basic_two_sequences(self):
        features = parse_tbl(_tbl(TBL_BASIC))
        assert set(features.keys()) == {"scaffold_1", "scaffold_2"}

    def test_feature_count_scaffold1(self):
        features = parse_tbl(_tbl(TBL_BASIC))
        # 2 coordinate lines + 2 qualifier lines
        assert len(features["scaffold_1"]) == 4

    def test_gene_coordinates(self):
        features = parse_tbl(_tbl(TBL_BASIC))
        first = features["scaffold_1"][0]
        assert first[0] == "1"
        assert first[1] == "1000"
        assert first[2] == "gene"

    def test_qualifier_stored(self):
        features = parse_tbl(_tbl(TBL_BASIC))
        quals = [f for f in features["scaffold_1"] if f[0] == ""]
        assert len(quals) >= 1
        assert "locus_tag\tGENE_001" in quals[0][3]

    def test_blank_lines_skipped(self):
        features = parse_tbl(_tbl(TBL_WITH_BLANKS))
        assert "scaffold_1" in features
        assert "scaffold_2" in features

    def test_comment_lines_skipped(self):
        features = parse_tbl(_tbl(TBL_WITH_COMMENT))
        assert "scaffold_1" in features

    def test_partial_coordinates(self):
        features = parse_tbl(_tbl(TBL_PARTIAL_COORDS))
        first = features["scaffold_1"][0]
        assert first[0] == "<1"
        assert first[1] == ">1000"


# ---------------------------------------------------------------------------
# parse_adjustments
# ---------------------------------------------------------------------------

ADJ_BASIC = """\
# Some header
#accession\toriginal_length\taction\tranges
scaffold_1\t1000\tACTION_TRIM\t1..50
scaffold_2\t2000\tACTION_TRIM\t1951..2000
"""

ADJ_INTERNAL = """\
#accession\toriginal_length\taction\tranges
scaffold_1\t1000\tACTION_TRIM\t200..400
"""

ADJ_NO_HEADER = """\
scaffold_1\t1000\tACTION_TRIM\t1..50
"""

ADJ_MULTI_RANGE = """\
#accession\toriginal_length\taction\tranges
scaffold_1\t1000\tACTION_TRIM\t1..30,971..1000
"""


class TestParseAdjustments:
    def test_basic_left_trim(self):
        adj = parse_adjustments(_adj(ADJ_BASIC))
        assert "scaffold_1" in adj
        # range 1..50 with start=1 → trim_left
        assert adj["scaffold_1"][0] == [1000, "ACTION_TRIM", 1, 50]

    def test_basic_right_trim(self):
        adj = parse_adjustments(_adj(ADJ_BASIC))
        assert "scaffold_2" in adj
        # range 1951..2000, end==original_length=2000 → trim_right
        assert adj["scaffold_2"][0] == [2000, "ACTION_TRIM", 1951, 2000]

    def test_skips_lines_before_header(self, capsys):
        adj = parse_adjustments(_adj(ADJ_NO_HEADER))
        # Without the #accession header the rows are skipped
        assert adj == {}

    def test_multi_range(self):
        adj = parse_adjustments(_adj(ADJ_MULTI_RANGE))
        # Two ranges for scaffold_1
        assert len(adj["scaffold_1"]) == 2


# ---------------------------------------------------------------------------
# fix_tbl (integration of parse_tbl + parse_adjustments + coordinate fixing)
# ---------------------------------------------------------------------------

TBL_FOR_FIX = """\
>Feature scaffold_1
1\t100\tgene
\t\t\tlocus_tag\tGENE_A
200\t300\tgene
\t\t\tlocus_tag\tGENE_B
"""

ADJ_LEFT_50 = """\
#accession\toriginal_length\taction\tranges
scaffold_1\t1000\tACTION_TRIM\t1..50
"""

ADJ_RIGHT_800 = """\
#accession\toriginal_length\taction\tranges
scaffold_1\t1000\tACTION_TRIM\t801..1000
"""

ADJ_NONE = """\
#accession\toriginal_length\taction\tranges
scaffold_2\t1000\tACTION_TRIM\t1..50
"""


class TestFixTbl:
    def test_trim_left_shifts_coordinates(self):
        out = _out()
        fix_tbl(_tbl(TBL_FOR_FIX), _adj(ADJ_LEFT_50), out)
        content = out.getvalue()
        lines = [l for l in content.splitlines() if "\t" in l and not l.startswith(">")]
        coord_lines = [l for l in lines if not l.startswith("\t")]
        # Original gene 1-100, after trim_left=50 → 1-50 (start clamped to 1)
        first_coords = coord_lines[0].split("\t")[:2]
        assert first_coords[0] == "1"    # max(1, 1-50)
        assert first_coords[1] == "50"   # 100-50

    def test_trim_left_second_gene_shifted(self):
        out = _out()
        fix_tbl(_tbl(TBL_FOR_FIX), _adj(ADJ_LEFT_50), out)
        content = out.getvalue()
        lines = [l for l in content.splitlines() if "\t" in l and not l.startswith(">")]
        coord_lines = [l for l in lines if not l.startswith("\t")]
        second_coords = coord_lines[2].split("\t")[:2]
        assert second_coords[0] == "150"  # 200-50
        assert second_coords[1] == "250"  # 300-50

    def test_trim_right_clamps_overlapping_feature(self):
        out = _out()
        fix_tbl(_tbl(TBL_FOR_FIX), _adj(ADJ_RIGHT_800), out)
        content = out.getvalue()
        lines = [l for l in content.splitlines() if "\t" in l and not l.startswith(">")]
        coord_lines = [l for l in lines if not l.startswith("\t")]
        # gene at 200-300: neither start(200) nor end(300) >= 801, so unchanged
        second_coords = coord_lines[2].split("\t")[:2]
        assert second_coords[0] == "200"
        assert second_coords[1] == "300"

    def test_no_adjustment_for_this_sequence(self):
        # ADJ_NONE applies to scaffold_2, not scaffold_1 in TBL_FOR_FIX
        out = _out()
        fix_tbl(_tbl(TBL_FOR_FIX), _adj(ADJ_NONE), out)
        content = out.getvalue()
        lines = [l for l in content.splitlines() if "\t" in l and not l.startswith(">")]
        coord_lines = [l for l in lines if not l.startswith("\t")]
        # Coordinates should be unchanged
        assert coord_lines[0].split("\t")[:2] == ["1", "100"]

    def test_output_has_feature_header(self):
        out = _out()
        fix_tbl(_tbl(TBL_FOR_FIX), _adj(ADJ_NONE), out)
        content = out.getvalue()
        assert ">Feature scaffold_1" in content

    def test_qualifiers_preserved(self):
        out = _out()
        fix_tbl(_tbl(TBL_FOR_FIX), _adj(ADJ_NONE), out)
        content = out.getvalue()
        assert "GENE_A" in content
        assert "GENE_B" in content
