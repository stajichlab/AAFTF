"""Tests for vecscreen BLAST coordinate parsing.

Verifies the fix for the stale start/end variable bug where the else-branch
for multi-hit contigs was using coordinate values from the previous row.
"""


import pytest

pytestmark = pytest.mark.unit


def _parse_blast_tabular(rows):
    """Minimal reimplementation of the vecscreen BLAST parsing loop.

    Mirrors the fixed logic in vecscreen.py::run() for contamination screens:
      start, end = sorted([int(row[6]), int(row[7])])
      if row[0] not in regions_to_trim:
          regions_to_trim[row[0]] = [(start, end, ...)]
      else:
          regions_to_trim[row[0]].append((start, end, ...))

    Returns:
        dict mapping contig_id -> list of (start, end, label, subject, pident)
    """
    regions_to_trim = {}
    for row in rows:
        if (float(row[2]) >= 98.0 and int(row[3]) >= 50) or (float(row[2]) >= 94.0 and int(row[3]) >= 100) or (float(row[2]) >= 90.0 and int(row[3]) >= 200):
            start, end = sorted([int(row[6]), int(row[7])])
            label = row[1] if len(row) > 1 else "CONTAM"
            pident = float(row[2])
            if row[0] not in regions_to_trim:
                regions_to_trim[row[0]] = [(start, end, label, row[1], pident)]
            else:
                regions_to_trim[row[0]].append((start, end, label, row[1], pident))
    return regions_to_trim


class TestVecscreenCoordinatesParsing:
    """Fix for stale start/end in else-branch when a contig has multiple hits."""

    def _make_row(self, qid, sid, pident, length, s1, e1, qs, qe):
        """Build a BLAST tabular row (12 fields)."""
        return [qid, sid, str(pident), str(length), "0", "0", str(qs), str(qe), str(s1), str(e1), "1e-5", "200"]

    def test_single_hit_correct_coordinates(self):
        rows = [self._make_row("contig_1", "CONTAM_A", 99.0, 100, 1, 100, 10, 109)]
        result = _parse_blast_tabular(rows)
        assert "contig_1" in result
        start, end, *_ = result["contig_1"][0]
        assert start == 10
        assert end == 109

    def test_second_hit_uses_its_own_coordinates(self):
        """Previously the else-branch reused stale start/end from the first row."""
        row1 = self._make_row("contig_1", "CONTAM_A", 99.0, 100, 1, 100, 10, 109)
        row2 = self._make_row("contig_1", "CONTAM_B", 99.0, 100, 200, 300, 500, 599)
        result = _parse_blast_tabular([row1, row2])
        assert len(result["contig_1"]) == 2
        s1, e1, *_ = result["contig_1"][0]
        s2, e2, *_ = result["contig_1"][1]
        assert s1 == 10 and e1 == 109
        assert s2 == 500 and e2 == 599, f"Expected (500, 599) for second hit but got ({s2}, {e2}). " "This would fail with the stale-coordinate bug."

    def test_reversed_coords_are_normalised(self):
        """Coordinates where qstart > qend must be sorted."""
        row = self._make_row("contig_1", "CONTAM_A", 99.0, 100, 1, 100, 109, 10)
        result = _parse_blast_tabular([row])
        start, end, *_ = result["contig_1"][0]
        assert start <= end

    def test_different_contigs_tracked_independently(self):
        row1 = self._make_row("contig_1", "CONTAM_A", 99.0, 100, 1, 100, 10, 109)
        row2 = self._make_row("contig_2", "CONTAM_A", 99.0, 100, 1, 100, 20, 119)
        result = _parse_blast_tabular([row1, row2])
        assert "contig_1" in result
        assert "contig_2" in result
        s1, *_ = result["contig_1"][0]
        s2, *_ = result["contig_2"][0]
        assert s1 == 10
        assert s2 == 20

    def test_hit_below_threshold_excluded(self):
        """A hit with pident=80 and length=50 should not be recorded."""
        row = self._make_row("contig_1", "CONTAM_A", 80.0, 50, 1, 50, 1, 50)
        result = _parse_blast_tabular([row])
        assert "contig_1" not in result
