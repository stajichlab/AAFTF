"""Unit tests for pure-Python functions in AAFTF/depth.py.

Tests cover FASTQ counting, mosdepth output parsing, and the outlier
detection arithmetic.  No external tools (minimap2, samtools, mosdepth)
are required.
"""

import gzip
import math
import os

import pytest

from AAFTF.depth import (
    _coverage_breadth_from_dist,
    count_fastq_reads,
    parse_mosdepth_summary,
)

from tests.conftest import MOSDEPTH_DIST, MOSDEPTH_SUMMARY, make_fastq_text

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# count_fastq_reads
# ---------------------------------------------------------------------------

class TestCountFastqReads:
    def test_plain_fastq(self, fastq_file):
        assert count_fastq_reads(str(fastq_file)) == 10

    def test_gz_fastq(self, gz_fastq_file):
        assert count_fastq_reads(str(gz_fastq_file)) == 10

    def test_single_read(self, tmp_path):
        p = tmp_path / "one.fastq"
        p.write_text(make_fastq_text(1))
        assert count_fastq_reads(str(p)) == 1

    def test_50_reads(self, tmp_path):
        p = tmp_path / "fifty.fastq"
        p.write_text(make_fastq_text(50))
        assert count_fastq_reads(str(p)) == 50

    def test_gz_50_reads(self, tmp_path):
        p = tmp_path / "fifty.fastq.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(make_fastq_text(50))
        assert count_fastq_reads(str(p)) == 50

    def test_returns_minus_one_on_error(self, tmp_path):
        # Corrupted gzip content → error → returns -1
        p = tmp_path / "bad.fastq.gz"
        p.write_bytes(b"not a valid gzip file")
        result = count_fastq_reads(str(p))
        assert result == -1


# ---------------------------------------------------------------------------
# parse_mosdepth_summary
# ---------------------------------------------------------------------------

class TestParseMosdepthSummary:
    def test_contig_count(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        # 9 nuclear + 1 outlier = 10 contig rows
        assert len(contigs) == 10

    def test_total_row_present(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        assert total is not None

    def test_total_mean(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        assert abs(total["mean"] - 20.88) < 0.01

    def test_total_length(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        # 9 × 10000 + 1000 = 91000
        assert total["length"] == 91_000

    def test_nuclear_contig_mean(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        nuclear = [c for c in contigs if c["chrom"].startswith("scaffold_")
                   and c["chrom"] != "scaffold_outlier"]
        assert all(abs(c["mean"] - 10.0) < 0.01 for c in nuclear)

    def test_outlier_contig_mean(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        outlier = next(c for c in contigs if c["chrom"] == "scaffold_outlier")
        assert abs(outlier["mean"] - 1000.0) < 0.01

    def test_no_total_row_excluded_from_contigs(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        assert all(c["chrom"] != "total" for c in contigs)

    def test_handles_missing_total_row(self, tmp_path):
        # Summary with no 'total' line
        content = (
            "chrom\tlength\tbases\tmean\tmin\tmax\n"
            "scaffold_1\t1000\t50000\t50.0\t0\t100\n"
        )
        p = tmp_path / "no_total.txt"
        p.write_text(content)
        total, contigs = parse_mosdepth_summary(str(p))
        assert total is None
        assert len(contigs) == 1


# ---------------------------------------------------------------------------
# _coverage_breadth_from_dist
# ---------------------------------------------------------------------------

class TestCoverageBreadthFromDist:
    def test_returns_correct_fraction(self, mosdepth_dist_file):
        result = _coverage_breadth_from_dist(str(mosdepth_dist_file))
        assert result is not None
        assert abs(result - 0.98) < 1e-6

    def test_returns_none_for_missing_file(self, tmp_path):
        result = _coverage_breadth_from_dist(str(tmp_path / "nonexistent.txt"))
        assert result is None

    def test_returns_none_when_threshold_1_absent(self, tmp_path):
        # File has threshold 0 and 2 but not 1
        content = "total\t0\t1.00000\ntotal\t2\t0.96000\n"
        p = tmp_path / "no_thresh1.txt"
        p.write_text(content)
        result = _coverage_breadth_from_dist(str(p))
        assert result is None


# ---------------------------------------------------------------------------
# Outlier detection arithmetic (mirrors logic in depth.run())
# ---------------------------------------------------------------------------

class TestOutlierDetectionMath:
    """Verify the SD-based flagging formula used in depth.run()."""

    def _compute(self, depths, mean_depth):
        sd = math.sqrt(
            sum((d - mean_depth) ** 2 for d in depths) / len(depths)
        )
        threshold = mean_depth + 3.0 * sd
        outliers = [d for d in depths if d > threshold]
        return sd, threshold, outliers

    def test_outlier_flagged(self, mosdepth_summary_file):
        # Use the test fixture data: 9 × 10x, 1 × 1000x, mean=20.88
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        depths = [c["mean"] for c in contigs]
        mean_depth = total["mean"]   # 20.88 (length-weighted)
        sd, threshold, outliers = self._compute(depths, mean_depth)
        # scaffold_outlier (1000x) must be above the threshold
        assert 1000.0 in outliers

    def test_no_outlier_when_uniform(self):
        depths = [50.0] * 10
        sd, threshold, outliers = self._compute(depths, 50.0)
        assert sd == 0.0
        assert outliers == []

    def test_threshold_formula(self):
        # Simple: mean=10, all depths=10 except one at 1000
        depths = [10.0] * 9 + [1000.0]
        mean = 10.0
        sd, threshold, outliers = self._compute(depths, mean)
        expected_sd = math.sqrt((9 * 0 + 990 ** 2) / 10)
        assert abs(sd - expected_sd) < 0.01
        assert threshold == pytest.approx(mean + 3 * expected_sd)
        assert 1000.0 in outliers
