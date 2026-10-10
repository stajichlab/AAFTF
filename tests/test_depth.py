"""Unit tests for pure-Python functions in aaftf/depth.py.

Tests cover FASTQ counting, mosdepth output parsing, outlier detection
arithmetic, quantize bin parsing, plot prefix derivation, and quantized
BED reading.  No external tools (minimap2, samtools, mosdepth) are required.
"""

import gzip
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.depth import (
    _coverage_breadth_from_dist,
    _depth_outliers,
    _get_plot_prefix,
    _paginate_by_length_ratio,
    _parse_quantize_bins,
    _read_quantized_bed,
    map_reads,
    parse_mosdepth_summary,
    run,
    run_flagstat,
)

pytestmark = pytest.mark.unit


class TestRunGuards:
    def test_longreads_without_preset_raises(self, tmp_path):
        asm = tmp_path / "asm.fa"
        asm.write_text(">contig1\nACGT\n")
        lr = tmp_path / "lr.fq"
        lr.write_text("@r\nACGT\n+\nIIII\n")
        with pytest.raises(ValueError):
            run(input=str(asm), longreads=str(lr))


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
        content = "chrom\tlength\tbases\tmean\tmin\tmax\nscaffold_1\t1000\t50000\t50.0\t0\t100\n"
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
# _depth_outliers (statistics used by depth.run())
# ---------------------------------------------------------------------------


def _row(name, length, mean):
    return {"chrom": name, "length": length, "mean": mean}


class TestDepthOutliers:
    def test_outlier_flagged(self, mosdepth_summary_file):
        total, contigs = parse_mosdepth_summary(str(mosdepth_summary_file))
        stats = _depth_outliers(contigs, total, 500)
        assert stats["mean_depth"] == pytest.approx(20.88)
        assert stats["n_outliers_3sd"] == 1
        assert stats["threshold_3sd"] < 1000.0

    def test_uniform_depths_nothing_flagged(self):
        rows = [_row(f"c{i}", 1000, 50.0) for i in range(10)]
        stats = _depth_outliers(rows, None, 500)
        assert stats["sd_depth"] == 0.0
        assert stats["threshold_3sd"] == 50.0
        assert stats["n_outliers_2sd"] == 0
        assert stats["n_outliers_3sd"] == 0

    def test_population_sd(self):
        rows = [_row(f"c{i}", 1000, 10.0) for i in range(9)] + [_row("big", 1000, 1000.0)]
        stats = _depth_outliers(rows, {"length": 10000, "mean": 10.0}, 500)
        # population SD (divide by n=10) around the mosdepth mean of 10
        assert stats["sd_depth"] == pytest.approx(990 / 10**0.5)
        assert stats["contig_arith_mean"] == pytest.approx(109.0)

    def test_short_contig_excluded_from_stats(self):
        rows = [_row("a", 1000, 10.0), _row("b", 1000, 20.0), _row("tiny", 100, 5000.0)]
        stats = _depth_outliers(rows, None, 500)
        assert stats["contig_arith_mean"] == pytest.approx(15.0)
        assert stats["sd_depth"] == pytest.approx(5.0)
        # the short contig is still counted as an outlier
        assert stats["n_outliers_3sd"] == 1

    def test_no_total_row_uses_contig_mean(self):
        rows = [_row("a", 1000, 10.0), _row("b", 1000, 30.0)]
        stats = _depth_outliers(rows, None, 500)
        assert stats["mosdepth_mean_depth"] is None
        assert stats["mean_depth"] == pytest.approx(20.0)

    def test_total_row_preferred_over_contig_mean(self):
        rows = [_row("a", 1000, 10.0), _row("b", 1000, 30.0)]
        stats = _depth_outliers(rows, {"length": 2000, "mean": 12.0}, 500)
        assert stats["mean_depth"] == 12.0
        assert stats["contig_arith_mean"] == pytest.approx(20.0)

    def test_no_contigs(self):
        stats = _depth_outliers([], None, 500)
        assert stats == {
            "mosdepth_mean_depth": None,
            "contig_arith_mean": 0.0,
            "mean_depth": 0.0,
            "sd_depth": 0.0,
            "threshold_2sd": 0.0,
            "threshold_3sd": 0.0,
            "n_outliers_2sd": 0,
            "n_outliers_3sd": 0,
        }

    def test_all_contigs_short_falls_back_to_total(self):
        stats = _depth_outliers([_row("a", 100, 7.0)], {"length": 100, "mean": 7.0}, 500)
        assert stats["mean_depth"] == 7.0
        assert stats["sd_depth"] == 0.0
        # thresholds are 0, so the contig counts as a 3SD outlier
        assert stats["n_outliers_3sd"] == 1

    def test_2sd_vs_3sd_counts(self):
        # 8 contigs at 0 and 2 at 10: mean 2, population SD 4 -> 2SD=10, 3SD=14
        rows = [_row(f"z{i}", 1000, 0.0) for i in range(8)] + [_row(f"t{i}", 1000, 10.0) for i in range(2)]
        stats = _depth_outliers(rows, None, 500)
        assert stats["threshold_2sd"] == pytest.approx(10.0)
        assert stats["threshold_3sd"] == pytest.approx(14.0)
        assert stats["n_outliers_2sd"] == 0  # 10 is not > 10
        # add contigs in (10, 14] and > 14
        rows += [_row("mid", 1000, 12.0), _row("high", 1000, 100.0)]
        stats = _depth_outliers(rows, {"length": 12000, "mean": 2.0}, 500)
        t2, t3 = stats["threshold_2sd"], stats["threshold_3sd"]
        assert stats["n_outliers_2sd"] == sum(1 for r in rows if t2 < r["mean"] <= t3)
        assert stats["n_outliers_3sd"] == sum(1 for r in rows if r["mean"] > t3) >= 1


class TestMapReadsErrors:
    def _call(self, tmp_path, **kw):
        args = dict(
            genome=str(tmp_path / "asm.fa"),
            read1="r1.fq",
            read2=None,
            longreads="lr.fq",
            workdir=str(tmp_path),
            cpus=1,
            illumina_preset="sr",
            longread_preset="map-ont",
            aligner="minimap2",
            debug=False,
        )
        args.update(kw)
        return map_reads(**args)

    def test_merge_failure_raises(self, tmp_path):
        with patch("aaftf.depth.align_to_sorted_bam"), patch("aaftf.depth.run_cmd", side_effect=lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, 1 if cmd[1] == "merge" else 0)):
            with pytest.raises(RuntimeError, match="merge"):
                self._call(tmp_path)

    def test_index_failure_raises(self, tmp_path):
        with patch("aaftf.depth.align_to_sorted_bam"), patch("aaftf.depth.run_cmd", side_effect=lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, 1 if cmd[1] == "index" else 0)):
            with pytest.raises(RuntimeError):
                self._call(tmp_path)

    def test_merge_success_returns_combined(self, tmp_path):
        with patch("aaftf.depth.align_to_sorted_bam"), patch("aaftf.depth.run_cmd", side_effect=lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, 0)):
            _, _, combined = self._call(tmp_path)
        assert combined.endswith("combined.sorted.bam")

    def test_longreads_without_preset_raises(self, tmp_path):
        with patch("aaftf.depth.align_to_sorted_bam"):
            with pytest.raises(ValueError):
                self._call(tmp_path, read1=None, longread_preset=None)


class TestRunFlagstat:
    def test_nonzero_raises(self):
        fake = subprocess.CompletedProcess(["samtools"], 1, stdout="", stderr="boom")
        with patch("aaftf.depth.subprocess.run", return_value=fake):
            with pytest.raises(RuntimeError, match="boom"):
                run_flagstat("x.bam")

    def test_success_returns_stdout(self):
        fake = subprocess.CompletedProcess(["samtools"], 0, stdout="10 + 0 mapped\n", stderr="")
        with patch("aaftf.depth.subprocess.run", return_value=fake):
            assert run_flagstat("x.bam") == "10 + 0 mapped\n"


# ---------------------------------------------------------------------------
# _parse_quantize_bins
# ---------------------------------------------------------------------------


class TestParseQuantizeBins:
    def test_default_returns_five_labels(self):
        labels, colors = _parse_quantize_bins()
        assert len(labels) == 5

    def test_default_uses_canonical_names(self):
        labels, colors = _parse_quantize_bins()
        assert labels[0] == "NO_COVERAGE"
        assert labels[2] == "CALLABLE"
        assert labels[4] == "VERY_HIGH_COVERAGE"

    def test_default_labels_get_canonical_colors(self):
        labels, colors = _parse_quantize_bins()
        assert colors["NO_COVERAGE"] == "#000080"
        assert colors["VERY_HIGH_COVERAGE"] == "#800000"


# ---------------------------------------------------------------------------
# _get_plot_prefix
# ---------------------------------------------------------------------------


class TestGetPlotPrefix:
    def test_strips_fasta_extension(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("genome.fasta", report)
        assert Path(prefix).name == "genome"

    def test_strips_fa_extension(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("assembly.fa", report)
        assert Path(prefix).name == "assembly"

    def test_strips_fasta_gz_extension(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("genome.fasta.gz", report)
        assert Path(prefix).name == "genome"

    def test_prefix_dir_matches_report_dir(self, tmp_path):
        report = str(tmp_path / "results" / "coverage_stats.txt")
        Path(tmp_path / "results").mkdir(parents=True, exist_ok=True)
        prefix = _get_plot_prefix("genome.fasta", report)
        assert Path(prefix).parent == tmp_path / "results"

    def test_dotted_name_preserves_stem(self, tmp_path):
        report = str(tmp_path / "coverage_stats.txt")
        prefix = _get_plot_prefix("strain.final.sorted.fasta", report)
        assert Path(prefix).name == "strain.final.sorted"


# ---------------------------------------------------------------------------
# _paginate_by_length_ratio
# ---------------------------------------------------------------------------


class TestPaginateByLengthRatio:
    def _rows(self, lengths):
        return [{"chrom": f"ctg{i}", "length": length} for i, length in enumerate(lengths)]

    def test_similar_lengths_stay_on_one_page(self):
        pages = _paginate_by_length_ratio(self._rows([1000, 900, 800, 700]))
        assert len(pages) == 1
        assert pages[0] == ["ctg0", "ctg1", "ctg2", "ctg3"]

    def test_splits_when_ratio_exceeded(self):
        # 1000 / 50 = 20 > default max_ratio (10), so ctg2 starts a new page
        pages = _paginate_by_length_ratio(self._rows([1000, 200, 50]))
        assert pages == [["ctg0", "ctg1"], ["ctg2"]]

    def test_splits_when_max_per_page_exceeded(self):
        pages = _paginate_by_length_ratio(self._rows([100, 100, 100]), max_per_page=2)
        assert pages == [["ctg0", "ctg1"], ["ctg2"]]

    def test_empty_input_returns_no_pages(self):
        assert _paginate_by_length_ratio([]) == []


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
