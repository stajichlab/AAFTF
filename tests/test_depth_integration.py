"""Integration tests for AAFTF/depth.py.

These tests invoke the *real* external tools (minimap2, bwa, samtools and
mosdepth) on a small synthetic genome with reads sampled directly from it.
They exercise the actual mapping + coverage pipeline end-to-end rather than
the pure-Python parsing covered by tests/test_depth.py.

The synthetic genome contains four "nuclear" contigs sampled to ~30x and one
short "organelle" contig sampled to ~400x, so outlier detection has a real
high-depth contig to find.

Run with:
    pixi run pytest tests/test_depth_integration.py -m integration -v

Tests are skipped automatically when the required tools are not on PATH, so
the suite still passes in a bare environment.
"""

import gzip
import random
import shutil
from argparse import Namespace
from pathlib import Path

import pytest

from AAFTF import depth

pytestmark = pytest.mark.integration

# ---------------------------------------------------------------------------
# Tool availability — skip the whole module if the core tools are missing
# ---------------------------------------------------------------------------

_have_minimap2 = shutil.which("minimap2") is not None
_have_samtools = shutil.which("samtools") is not None
_have_mosdepth = shutil.which("mosdepth") is not None
_have_bwa = shutil.which("bwa") is not None

_missing = [
    name
    for name, present in (
        ("minimap2", _have_minimap2),
        ("samtools", _have_samtools),
        ("mosdepth", _have_mosdepth),
    )
    if not present
]

requires_tools = pytest.mark.skipif(
    bool(_missing),
    reason="missing required tools: " + ", ".join(_missing),
)
requires_bwa = pytest.mark.skipif(not _have_bwa, reason="bwa not in PATH")

# ---------------------------------------------------------------------------
# Synthetic genome / reads
# ---------------------------------------------------------------------------

# Contigs: (name, length_bp, target_mean_depth)
# Many nuclear contigs at ~30x keep the population SD small enough that the
# short, high-depth organelle (400x) lands above the mean + 3*SD threshold and
# is reliably flagged as an outlier.
_CONTIG_SPEC = [(f"chr{i}", 3000, 30) for i in range(1, 21)]
_CONTIG_SPEC.append(("organelle", 1000, 400))
_READ_LEN = 100
_FRAG_LEN = 300

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _random_seq(rng, length):
    """Low-complexity-avoiding random ACGT sequence for unique mapping."""
    return "".join(rng.choice("ACGT") for _ in range(length))


def _build_genome_and_reads(out_dir):
    """Create a synthetic genome plus paired reads sampled from it.

    Returns (genome_path, r1_path, r2_path).  Reads are exact substrings of
    the genome (FR orientation) so they map near-perfectly with minimap2 -ax sr.
    """
    rng = random.Random(1234)

    genome = {}
    for name, length, _depth in _CONTIG_SPEC:
        genome[name] = _random_seq(rng, length)

    genome_path = out_dir / "genome.fasta"
    with open(genome_path, "w") as fh:
        for name, _length, _depth in _CONTIG_SPEC:
            seq = genome[name]
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")

    # Build read pairs from fragments to hit the target per-contig depth.
    r1_records = []
    r2_records = []
    read_idx = 0
    for name, length, target_depth in _CONTIG_SPEC:
        seq = genome[name]
        n_frags = max(1, target_depth * length // (2 * _READ_LEN))
        max_start = length - _FRAG_LEN
        for _ in range(n_frags):
            start = rng.randint(0, max_start)
            frag = seq[start : start + _FRAG_LEN]
            fwd = frag[:_READ_LEN]
            rev = _revcomp(frag[-_READ_LEN:])
            qual = "I" * _READ_LEN
            r1_records.append(f"@read{read_idx}/1\n{fwd}\n+\n{qual}\n")
            r2_records.append(f"@read{read_idx}/2\n{rev}\n+\n{qual}\n")
            read_idx += 1

    r1_path = out_dir / "reads_R1.fq.gz"
    r2_path = out_dir / "reads_R2.fq.gz"
    with gzip.open(r1_path, "wt") as fh:
        fh.write("".join(r1_records))
    with gzip.open(r2_path, "wt") as fh:
        fh.write("".join(r2_records))

    return genome_path, r1_path, r2_path


@pytest.fixture(scope="module")
def synthetic_data(tmp_path_factory):
    """Module-scoped synthetic genome + reads (built once, reused by tests)."""
    out_dir = tmp_path_factory.mktemp("depth_synth")
    genome, r1, r2 = _build_genome_and_reads(out_dir)
    return Namespace(genome=str(genome), r1=str(r1), r2=str(r2))


def _depth_args(synthetic_data, tmp_path, **overrides):
    """Build a complete Namespace for depth.run() with sensible defaults."""
    defaults = dict(
        input=synthetic_data.genome,
        left=synthetic_data.r1,
        right=synthetic_data.r2,
        longreads=None,
        illumina_preset="sr",
        longread_preset="map-ont",
        aligner="minimap2",
        cpus=2,
        workdir=str(tmp_path / "workdir"),
        out=str(tmp_path / "coverage_stats.txt"),
        min_contig_len=500,
        no_plot=True,
        plot_format="pdf",
        quantize=depth._DEFAULT_QUANTIZE,
        quantize_labels=None,
        debug=False,
        pipe=True,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


# ---------------------------------------------------------------------------
# Low-level pipeline functions against real tool output
# ---------------------------------------------------------------------------


@requires_tools
class TestMapReadsReal:
    def test_produces_sorted_indexed_bam(self, synthetic_data, tmp_path):
        workdir = tmp_path / "wd"
        workdir.mkdir()
        bam_ill, bam_lr, bam_combined = depth.map_reads(
            genome=synthetic_data.genome,
            reads_left=synthetic_data.r1,
            reads_right=synthetic_data.r2,
            longreads=None,
            workdir=str(workdir),
            cpus=2,
            illumina_preset="sr",
            longread_preset="map-ont",
            aligner="minimap2",
            debug=False,
        )
        assert bam_lr is None
        assert bam_combined == bam_ill
        assert Path(bam_combined).exists()
        assert Path(bam_combined).stat().st_size > 0
        assert Path(bam_combined + ".bai").exists()

    def test_flagstat_reports_high_mapping_rate(self, synthetic_data, tmp_path):
        workdir = tmp_path / "wd"
        workdir.mkdir()
        _, _, bam = depth.map_reads(
            genome=synthetic_data.genome,
            reads_left=synthetic_data.r1,
            reads_right=synthetic_data.r2,
            longreads=None,
            workdir=str(workdir),
            cpus=2,
            illumina_preset="sr",
            longread_preset="map-ont",
            aligner="minimap2",
            debug=False,
        )
        flagstat = depth.run_flagstat(bam)
        assert "mapped" in flagstat
        # reads are exact genome substrings: mapping rate should be near 100%
        line = next(ln for ln in flagstat.splitlines() if "mapped (" in ln)
        pct = float(line.split("(")[1].split("%")[0])
        assert pct > 95.0, f"mapping rate unexpectedly low: {line}"


@requires_tools
class TestMosdepthReal:
    def test_summary_lists_all_contigs(self, synthetic_data, tmp_path):
        workdir = tmp_path / "wd"
        workdir.mkdir()
        _, _, bam = depth.map_reads(
            genome=synthetic_data.genome,
            reads_left=synthetic_data.r1,
            reads_right=synthetic_data.r2,
            longreads=None,
            workdir=str(workdir),
            cpus=2,
            illumina_preset="sr",
            longread_preset="map-ont",
            aligner="minimap2",
            debug=False,
        )
        summary = depth.run_mosdepth(bam, str(workdir), cpus=2)
        assert Path(summary).exists()
        total, contigs = depth.parse_mosdepth_summary(summary)
        names = {c["chrom"] for c in contigs}
        assert {"chr1", "chr2", "chr3", "organelle"}.issubset(names)
        assert total is not None
        assert total["mean"] > 0

    def test_organelle_has_higher_depth(self, synthetic_data, tmp_path):
        workdir = tmp_path / "wd"
        workdir.mkdir()
        _, _, bam = depth.map_reads(
            genome=synthetic_data.genome,
            reads_left=synthetic_data.r1,
            reads_right=synthetic_data.r2,
            longreads=None,
            workdir=str(workdir),
            cpus=2,
            illumina_preset="sr",
            longread_preset="map-ont",
            aligner="minimap2",
            debug=False,
        )
        summary = depth.run_mosdepth(bam, str(workdir), cpus=2)
        _, contigs = depth.parse_mosdepth_summary(summary)
        by_name = {c["chrom"]: c["mean"] for c in contigs}
        nuclear_mean = sum(by_name[c] for c in ("chr1", "chr2", "chr3")) / 3
        # organelle was sampled to ~400x vs ~30x nuclear → clear separation
        assert by_name["organelle"] > 3 * nuclear_mean

    def test_quantized_mode_writes_bed(self, synthetic_data, tmp_path):
        workdir = tmp_path / "wd"
        workdir.mkdir()
        _, _, bam = depth.map_reads(
            genome=synthetic_data.genome,
            reads_left=synthetic_data.r1,
            reads_right=synthetic_data.r2,
            longreads=None,
            workdir=str(workdir),
            cpus=2,
            illumina_preset="sr",
            longread_preset="map-ont",
            aligner="minimap2",
            debug=False,
        )
        labels, env_dict = depth._parse_quantize_bins(depth._DEFAULT_QUANTIZE)
        summary, bed = depth.run_mosdepth_quantized(bam, str(workdir), cpus=2, quantize_str=depth._DEFAULT_QUANTIZE, env_dict=env_dict)
        assert Path(summary).exists()
        assert Path(bed).exists()
        contig_data = depth._read_quantized_bed(bed)
        assert "chr1" in contig_data
        # every interval label must come from our requested label set
        seen = {lbl for ivs in contig_data.values() for _, _, lbl in ivs}
        assert seen.issubset(set(labels))


# ---------------------------------------------------------------------------
# Full depth.run() end-to-end
# ---------------------------------------------------------------------------


@requires_tools
class TestDepthRunEndToEnd:
    def test_report_written_with_expected_sections(self, synthetic_data, tmp_path):
        args = _depth_args(synthetic_data, tmp_path)
        depth.run(None, args)

        report = Path(args.out)
        assert report.exists() and report.stat().st_size > 0
        text = report.read_text()
        assert "AAFTF Coverage Statistics Report" in text
        assert "Whole-Assembly Coverage" in text
        assert "Per-Contig Coverage" in text
        for contig in ("chr1", "chr2", "chr3", "organelle"):
            assert contig in text

    def test_organelle_flagged_as_outlier(self, synthetic_data, tmp_path):
        args = _depth_args(synthetic_data, tmp_path)
        depth.run(None, args)
        text = Path(args.out).read_text()
        # The organelle contig line should carry the OUTLIER flag.
        organelle_line = next(ln for ln in text.splitlines() if ln.strip().startswith("organelle"))
        assert "OUTLIER" in organelle_line, organelle_line

    def test_coverage_breadth_near_complete(self, synthetic_data, tmp_path):
        args = _depth_args(synthetic_data, tmp_path)
        depth.run(None, args)
        text = Path(args.out).read_text()
        breadth_line = next((ln for ln in text.splitlines() if "Bases covered" in ln), None)
        assert breadth_line is not None
        pct = float(breadth_line.split(":")[1].strip().rstrip("%"))
        assert pct > 95.0, breadth_line

    def test_workdir_removed_when_not_debug(self, synthetic_data, tmp_path, monkeypatch):
        # No explicit workdir → run() creates a uuid temp dir in CWD and is
        # expected to remove it when not in debug mode.  chdir into tmp_path so
        # the auto-created dir never lands in the repo.
        monkeypatch.chdir(tmp_path)
        args = _depth_args(synthetic_data, tmp_path, workdir=None)
        depth.run(None, args)
        assert Path(args.out).exists()
        leaked = list(tmp_path.glob("aaftf-depth_*"))
        assert not leaked, f"workdir not cleaned up: {leaked}"

    @pytest.mark.skipif(not depth.HAS_MATPLOTLIB, reason="matplotlib not installed")
    def test_plots_generated_when_enabled(self, synthetic_data, tmp_path):
        args = _depth_args(synthetic_data, tmp_path, no_plot=False, plot_format="png")
        depth.run(None, args)
        prefix = depth._get_plot_prefix(args.input, args.out)
        assert Path(prefix + ".depth_heatmap.png").exists()
        assert Path(prefix + ".depth_barplot.png").exists()
        assert Path(prefix + ".depth_histogram.png").exists()


@requires_tools
@requires_bwa
class TestDepthRunBwaAligner:
    def test_bwa_aligner_produces_report(self, synthetic_data, tmp_path):
        args = _depth_args(synthetic_data, tmp_path, aligner="bwa")
        depth.run(None, args)
        report = Path(args.out)
        assert report.exists() and report.stat().st_size > 0
        text = report.read_text()
        assert "organelle" in text
