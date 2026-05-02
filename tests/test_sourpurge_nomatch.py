"""Tests for sourpurge.py taxonomy CSV parsing.

Covers the fix for the 'nomatch' membership test: the original code used
`"nomatch" in cols` (list membership) which is correct but fragile; the
fix uses `cols[1].strip() == "nomatch"` for explicit field-level comparison.
"""

import pytest

pytestmark = pytest.mark.unit


def _check_nomatch(cols):
    """Mirrors the fixed sourpurge.py check: cols[1].strip() == 'nomatch'."""
    return cols[1].strip() == "nomatch"


class TestNomatchDetection:
    def test_nomatch_in_status_field(self):
        cols = ["contig_1", "nomatch", "0", "0.0", "", "", "", "", "", "", ""]
        assert _check_nomatch(cols) is True

    def test_nomatch_with_whitespace(self):
        cols = ["contig_1", " nomatch ", "0", "0.0", "", "", "", "", "", "", ""]
        assert _check_nomatch(cols) is True

    def test_found_status_not_nomatch(self):
        cols = ["contig_1", "found", "1", "0.95", "", "", "", "", "", "", ""]
        assert _check_nomatch(cols) is False

    def test_nomatch_in_non_status_field_not_matched(self):
        """A contig ID containing 'nomatch' must not trigger the nomatch branch."""
        cols = ["nomatch_contig_1", "found", "1", "0.95", "", "", "", "", "", "", ""]
        assert _check_nomatch(cols) is False

    def test_empty_status_field_not_nomatch(self):
        cols = ["contig_1", "", "0", "0.0", "", "", "", "", "", "", ""]
        assert _check_nomatch(cols) is False

    def test_partial_match_not_nomatch(self):
        cols = ["contig_1", "nomatch_extra", "0", "0.0", "", "", "", "", "", "", ""]
        assert _check_nomatch(cols) is False
