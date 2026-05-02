"""Unit tests for pure-Python functions in AAFTF/depth.py.

Tests cover FASTQ counting, mosdepth output parsing, outlier detection
arithmetic, quantize bin parsing, plot prefix derivation, and quantized
BED reading.  No external tools (minimap2, samtools, mosdepth) are required.
"""

import gzip
import math
import os

import pytest

from AAFTF.depth import (
    _coverage_breadth_from_dist,
    _get_plot_prefix,
    _parse_quantize_bins,
    _read_quantized_bed,
    count_fastq_reads,
    parse_mosdepth_summary,
)
from tests.conftest import make_fastq_text

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
        nuclear = [c for c in contigs if c["chrom"].startswith("scaffold_") and c["chrom"] != "scaffold_outlier"]
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
        content = "chrom\tlength\tbases\tmean\tmin\tmax\n" "scaffold_1\t1000\t50000\t50.0\t0\t100\n"
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
        # Uses population SD — contigs are the full population, not a sample.
        sd = math.sqrt(sum((d - mean_depth) ** 2 for d in depths) / len(depths))
        threshold = mean_depth + 3.0 * sd
        outliers = [d for d in depths if d > threshold]
        return sd, threshold, outliers

    def test_outlier_flagged(self, mosdepth_summary_file):
        # Use the test fixture data: 9 × 10x, 1 × 1000x, mean=20.88
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        depths = [c["mean"] for c in contigs]
        mean_depth = total["mean"]  # 20.88 (length-weighted)
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
        expected_sd = math.sqrt((9 * 0 + 990**2) / 10)
        assert abs(sd - expected_sd) < 0.01
        assert threshold == pytest.approx(mean + 3 * expected_sd)
        assert 1000.0 in outliers


# ---------------------------------------------------------------------------
# _parse_quantize_bins
# ---------------------------------------------------------------------------


class TestParseQuantizeBins:
    def test_default_bins_returns_five_labels(self):
        labels, env_dict = _parse_quantize_bins("0:1:4:100:200:")
        assert len(labels) == 5

    def test_default_bins_uses_canonical_names(self):
        labels, _ = _parse_quantize_bins("0:1:4:100:200:")
        assert labels[0] == "NO_COVERAGE"
        assert labels[2] == "CALLABLE"
        assert labels[4] == "VERY_HIGH_COVERAGE"

    def test_env_dict_keys_are_mosdepth_q_vars(self):
        _, env_dict = _parse_quantize_bins("0:1:4:100:200:")
        assert set(env_dict.keys()) == {"MOSDEPTH_Q0", "MOSDEPTH_Q1", "MOSDEPTH_Q2", "MOSDEPTH_Q3", "MOSDEPTH_Q4"}

    def test_env_dict_values_match_labels(self):
        labels, env_dict = _parse_quantize_bins("0:1:4:100:200:")
        for i, lbl in enumerate(labels):
            assert env_dict[f"MOSDEPTH_Q{i}"] == lbl

    def test_custom_bins_auto_named(self):
        labels, _ = _parse_quantize_bins("0:5:50:")
        assert labels == ["BIN_0", "BIN_1", "BIN_2"]

    def test_custom_labels_override(self):
        labels, _ = _parse_quantize_bins("0:5:50:", labels_str="NONE,LOW,HIGH")
        assert labels == ["NONE", "LOW", "HIGH"]

    def test_too_few_custom_labels_get_padded(self):
        labels, _ = _parse_quantize_bins("0:1:4:100:200:", labels_str="A,B")
        assert len(labels) == 5
        assert labels[0] == "A"
        assert labels[1] == "B"
        assert labels[2] == "BIN_2"

    def test_too_many_custom_labels_get_truncated(self):
        labels, _ = _parse_quantize_bins("0:5:", labels_str="X,Y,Z,W")
        assert labels == ["X", "Y"]

    def test_trailing_colon_not_counted_as_bin(self):
        labels, _ = _parse_quantize_bins("0:1:")
        assert len(labels) == 2

    def test_three_bin_custom_string(self):
        labels, env_dict = _parse_quantize_bins("0:10:100:")
        assert len(labels) == 3
        assert "MOSDEPTH_Q2" in env_dict


# ---------------------------------------------------------------------------
# _get_plot_prefix
# ---------------------------------------------------------------------------


class TestGetPlotPrefix:
    def test_strips_fasta_extension(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("genome.fasta", report)
        assert os.path.basename(prefix) == "genome"

    def test_strips_fa_extension(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("assembly.fa", report)
        assert os.path.basename(prefix) == "assembly"

    def test_strips_fasta_gz_extension(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("genome.fasta.gz", report)
        assert os.path.basename(prefix) == "genome"

    def test_prefix_dir_matches_report_dir(self, tmp_path):
        report = str(tmp_path / "results" / "coverage_stats.txt")
        os.makedirs(str(tmp_path / "results"), exist_ok=True)
        prefix = _get_plot_prefix("genome.fasta", report)
        assert os.path.dirname(prefix) == str(tmp_path / "results")

    def test_dotted_name_preserves_stem(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("strain.final.sorted.fasta", report)
        assert os.path.basename(prefix) == "strain.final.sorted"


# ---------------------------------------------------------------------------
# _read_quantized_bed
# ---------------------------------------------------------------------------


class TestReadQuantizedBed:
    def _write_bed(self, path, lines, gz=True):
        content = "\n".join(lines) + "\n"
        if gz:
            with gzip.open(path, "wt") as fh:
                fh.write(content)
        else:
            with open(path, "w") as fh:
                fh.write(content)

    def test_reads_gzipped_bed(self, tmp_path):
        bed = str(tmp_path / "test.quantized.bed.gz")
        self._write_bed(
            bed,
            [
                "scaffold_1\t0\t5000\tNO_COVERAGE",
                "scaffold_1\t5000\t100000\tCALLABLE",
                "scaffold_2\t0\t200000\tCALLABLE",
            ],
        )
        data = _read_quantized_bed(bed)
        assert "scaffold_1" in data
        assert "scaffold_2" in data
        assert len(data["scaffold_1"]) == 2

    def test_reads_plain_bed(self, tmp_path):
        bed = str(tmp_path / "test.quantized.bed")
        self._write_bed(
            bed,
            [
                "scf1\t0\t1000\tNO_COVERAGE",
                "scf1\t1000\t5000\tCALLABLE",
            ],
            gz=False,
        )
        data = _read_quantized_bed(bed)
        assert len(data["scf1"]) == 2

    def test_interval_fields_parsed_correctly(self, tmp_path):
        bed = str(tmp_path / "test.quantized.bed.gz")
        self._write_bed(bed, ["scaffold_1\t100\t500\tHIGH_COVERAGE"])
        data = _read_quantized_bed(bed)
        start, end, label = data["scaffold_1"][0]
        assert start == 100
        assert end == 500
        assert label == "HIGH_COVERAGE"

    def test_multiple_contigs_separated(self, tmp_path):
        bed = str(tmp_path / "test.quantized.bed.gz")
        self._write_bed(
            bed,
            [
                "chr1\t0\t1000\tCALLABLE",
                "chr2\t0\t2000\tLOW_COVERAGE",
                "chr1\t1000\t3000\tHIGH_COVERAGE",
            ],
        )
        data = _read_quantized_bed(bed)
        assert len(data["chr1"]) == 2
        assert len(data["chr2"]) == 1

    def test_short_lines_skipped(self, tmp_path):
        bed = str(tmp_path / "test.quantized.bed.gz")
        self._write_bed(
            bed,
            [
                "scaffold_1\t0\t1000\tCALLABLE",
                "bad_line",
                "scaffold_1\t1000\t2000\tHIGH_COVERAGE",
            ],
        )
        data = _read_quantized_bed(bed)
        assert len(data["scaffold_1"]) == 2
