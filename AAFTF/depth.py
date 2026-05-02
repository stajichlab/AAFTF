"""Calculate depth of coverage of reads mapped to a genome assembly.

Maps Illumina paired-end and/or long reads to the genome assembly using
minimap2 (or bwa for Illumina reads), then runs mosdepth to compute
per-contig and whole-assembly coverage statistics.  Contigs with mean
depth more than 3 standard deviations above the assembly mean are
flagged as possible contaminants or organelles.

When plotting is enabled (default, requires matplotlib), mosdepth is run in
quantized mode and three additional plots are produced alongside the report:
  <prefix>.depth_heatmap.<format>   — per-scaffold coverage-class tiles
  <prefix>.depth_barplot.<format>   — stacked bar of coverage-class proportions
  <prefix>.depth_histogram.<format> — histogram + boxplot of per-scaffold depths
"""

import gzip
import math
import os
import shutil
import subprocess
import sys
import uuid

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Patch

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from AAFTF.utility import SafeRemove, checkfile, printCMD, status, which

# ---------------------------------------------------------------------------
# Constants for quantized coverage classes
# ---------------------------------------------------------------------------

_DEFAULT_QUANTIZE = "0:1:4:100:200:"
_DEFAULT_LABELS = [
    "NO_COVERAGE",
    "LOW_COVERAGE",
    "CALLABLE",
    "HIGH_COVERAGE",
    "VERY_HIGH_COVERAGE",
]
_COV_COLOURS = {
    "NO_COVERAGE": "#313695",
    "LOW_COVERAGE": "#74add1",
    "CALLABLE": "#ffffbf",
    "HIGH_COVERAGE": "#f46d43",
    "VERY_HIGH_COVERAGE": "#a50026",
}
# Fallback palette for non-standard bin labels
_FALLBACK_COLOURS = ["#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026", "#762a83", "#1b7837"]

_SCAFFOLDS_PER_PAGE = 50

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def count_fastq_reads(fastq_file):
    """Count the number of reads in a FASTQ file.

    Args:
        fastq_file: Path to FASTQ file (gzip-compressed or plain).

    Returns:
        Integer count of reads, or -1 on failure.
    """
    try:
        opener = gzip.open if fastq_file.endswith(".gz") else open
        n_lines = 0
        with opener(fastq_file, "rt") as fh:
            for _ in fh:
                n_lines += 1
        return n_lines // 4
    except Exception:
        return -1


def _open_devnull():
    """Return an open handle to /dev/null."""
    return open(os.devnull, "w")


def map_reads(genome, reads_left, reads_right, longreads, workdir, cpus, illumina_preset, longread_preset, aligner, debug):
    """Map reads to the genome assembly and produce sorted, indexed BAM files.

    Illumina reads are mapped with minimap2 (default) or bwa mem.
    Long reads are always mapped with minimap2.
    When both read types are provided, the two BAM files are merged.

    Args:
        genome: Absolute path to genome assembly FASTA.
        reads_left: Path to forward Illumina FASTQ (or None).
        reads_right: Path to reverse Illumina FASTQ (or None).
        longreads: Path to long-read FASTQ (or None).
        workdir: Directory for intermediate files.
        cpus: Number of CPU threads.
        illumina_preset: minimap2 preset for Illumina reads (e.g. 'sr').
        longread_preset: minimap2 preset for long reads (e.g. 'map-ont').
        aligner: Aligner for Illumina reads ('minimap2' or 'bwa').
        debug: If True, show stderr from sub-processes.

    Returns:
        Tuple (bam_illumina, bam_longreads, bam_combined).
        bam_illumina and bam_longreads are None when the respective read
        type was not provided.  bam_combined is the single BAM to use for
        mosdepth (merged when both types are present).
    """
    devnull = _open_devnull()
    stderr_dest = None if debug else devnull

    bam_illumina = None
    bam_longreads = None

    # --- Illumina reads ---
    if reads_left:
        bam_illumina = os.path.join(workdir, "illumina.sorted.bam")
        if aligner == "bwa":
            # Copy genome to workdir so BWA index files stay contained
            genome_local = os.path.join(workdir, os.path.basename(genome))
            shutil.copyfile(genome, genome_local)
            bwa_index_cmd = ["bwa", "index", genome_local]
            printCMD(bwa_index_cmd)
            subprocess.run(bwa_index_cmd, stderr=stderr_dest)
            map_cmd = ["bwa", "mem", "-t", str(cpus), genome_local, reads_left]
            if reads_right:
                map_cmd.append(reads_right)
        else:
            map_cmd = [
                "minimap2",
                "-ax",
                illumina_preset,
                "-t",
                str(cpus),
                genome,
                reads_left,
            ]
            if reads_right:
                map_cmd.append(reads_right)

        printCMD(map_cmd)
        p1 = subprocess.Popen(map_cmd, stdout=subprocess.PIPE, stderr=stderr_dest)
        sort_cmd = ["samtools", "sort", "-@", str(cpus), "-o", bam_illumina, "-"]
        printCMD(sort_cmd)
        p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stderr=stderr_dest)
        p1.stdout.close()
        p2.communicate()
        subprocess.run(["samtools", "index", bam_illumina], stderr=stderr_dest)

    # --- Long reads ---
    if longreads:
        bam_longreads = os.path.join(workdir, "longreads.sorted.bam")
        map_cmd = [
            "minimap2",
            "-ax",
            longread_preset,
            "-t",
            str(cpus),
            genome,
            longreads,
        ]
        printCMD(map_cmd)
        p1 = subprocess.Popen(map_cmd, stdout=subprocess.PIPE, stderr=stderr_dest)
        sort_cmd = ["samtools", "sort", "-@", str(cpus), "-o", bam_longreads, "-"]
        printCMD(sort_cmd)
        p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stderr=stderr_dest)
        p1.stdout.close()
        p2.communicate()
        subprocess.run(["samtools", "index", bam_longreads], stderr=stderr_dest)

    devnull.close()

    # --- Combine ---
    if bam_illumina and bam_longreads:
        bam_combined = os.path.join(workdir, "combined.sorted.bam")
        devnull2 = _open_devnull()
        merge_cmd = [
            "samtools",
            "merge",
            "-@",
            str(cpus),
            "-f",
            bam_combined,
            bam_illumina,
            bam_longreads,
        ]
        printCMD(merge_cmd)
        subprocess.run(merge_cmd, stderr=None if debug else devnull2)
        subprocess.run(["samtools", "index", bam_combined], stderr=None if debug else devnull2)
        devnull2.close()
    elif bam_illumina:
        bam_combined = bam_illumina
    else:
        bam_combined = bam_longreads

    return bam_illumina, bam_longreads, bam_combined


def run_flagstat(bam_file):
    """Run samtools flagstat on a BAM file.

    Args:
        bam_file: Path to indexed BAM file.

    Returns:
        String containing the flagstat output lines.
    """
    result = subprocess.run(
        ["samtools", "flagstat", bam_file],
        capture_output=True,
        text=True,
    )
    return result.stdout


def run_mosdepth(bam_file, workdir, cpus, prefix="coverage", debug=False):
    """Run mosdepth to calculate per-contig depth statistics.

    Args:
        bam_file: Absolute path to the sorted, indexed BAM file.
        workdir: Directory for mosdepth output files.
        cpus: Number of CPU threads.
        prefix: Filename prefix for mosdepth outputs.
        debug: If True, show mosdepth stderr.

    Returns:
        Path to the mosdepth summary text file.
    """
    devnull = _open_devnull()
    mosdepth_prefix = os.path.join(os.path.abspath(workdir), prefix)
    cmd = [
        "mosdepth",
        "--threads",
        str(cpus),
        "--no-abbrev",
        mosdepth_prefix,
        os.path.abspath(bam_file),
    ]
    printCMD(cmd)
    subprocess.run(cmd, stderr=None if debug else devnull)
    devnull.close()
    return mosdepth_prefix + ".mosdepth.summary.txt"


def parse_mosdepth_summary(summary_file):
    """Parse a mosdepth summary file into per-contig rows.

    The summary file produced by mosdepth has tab-separated columns:
    chrom, length, bases, mean, min, max.  The final row named 'total'
    gives assembly-wide aggregated statistics.

    Args:
        summary_file: Path to mosdepth summary text file.

    Returns:
        Tuple (total_row, contig_rows) where each element is a dict with
        keys 'chrom', 'length', 'bases', 'mean'.  total_row may be None.
    """
    contigs = []
    total = None
    with open(summary_file) as fh:
        next(fh)  # skip header line
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            row = {
                "chrom": parts[0],
                "length": int(parts[1]),
                "bases": int(parts[2]),
                "mean": float(parts[3]),
            }
            if parts[0] == "total":
                total = row
            else:
                contigs.append(row)
    return total, contigs


def _coverage_breadth_from_dist(dist_file):
    """Read coverage breadth at >= 1x from mosdepth global distribution file.

    Args:
        dist_file: Path to *.mosdepth.global.dist.txt.

    Returns:
        Float fraction (0–1) of bases with depth >= 1, or None if unavailable.
    """
    if not os.path.exists(dist_file):
        return None
    try:
        with open(dist_file) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[0] == "total" and parts[1] == "1":
                    return float(parts[2])
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Quantized mosdepth helpers
# ---------------------------------------------------------------------------


def _parse_quantize_bins(quantize_str, labels_str=None):
    """Parse a mosdepth quantize string and return bin labels plus env-var dict.

    Args:
        quantize_str: Colon-separated bin boundaries with trailing colon,
            e.g. "0:1:4:100:200:".
        labels_str: Optional comma-separated label names, one per bin.
            When None, default names are used for the default quantize string;
            otherwise bins are named BIN_0, BIN_1, etc.

    Returns:
        Tuple (labels, env_dict) where labels is a list of strings and
        env_dict maps MOSDEPTH_Q0..QN to the corresponding label strings.
    """
    thresholds = [t for t in quantize_str.split(":") if t.strip()]
    n_bins = len(thresholds)

    if labels_str:
        labels = [lbl.strip() for lbl in labels_str.split(",")]
    elif quantize_str == _DEFAULT_QUANTIZE:
        labels = _DEFAULT_LABELS[:]
    else:
        labels = [f"BIN_{i}" for i in range(n_bins)]

    # pad or truncate to match bin count
    while len(labels) < n_bins:
        labels.append(f"BIN_{len(labels)}")
    labels = labels[:n_bins]

    env_dict = {f"MOSDEPTH_Q{i}": labels[i] for i in range(n_bins)}
    return labels, env_dict


def run_mosdepth_quantized(bam_file, workdir, cpus, quantize_str, env_dict, prefix="coverage", debug=False):
    """Run mosdepth in quantized mode.

    Produces a summary file (same format as run_mosdepth) plus a
    quantized BED file showing which coverage class each genomic interval
    belongs to.

    Args:
        bam_file: Absolute path to the sorted, indexed BAM file.
        workdir: Directory for mosdepth output files.
        cpus: Number of CPU threads.
        quantize_str: Colon-separated quantize boundaries, e.g. "0:1:4:100:200:".
        env_dict: Dict of MOSDEPTH_Q? env vars to set (label names for BED).
        prefix: Filename prefix for mosdepth outputs.
        debug: If True, show mosdepth stderr.

    Returns:
        Tuple (summary_file, quantized_bed) — paths to the two output files.
    """
    devnull = _open_devnull()
    mosdepth_prefix = os.path.join(os.path.abspath(workdir), prefix)
    run_env = os.environ.copy()
    run_env.update(env_dict)
    cmd = [
        "mosdepth",
        "--threads",
        str(cpus),
        "--no-abbrev",
        "-n",
        "--quantize",
        quantize_str,
        mosdepth_prefix,
        os.path.abspath(bam_file),
    ]
    printCMD(cmd)
    subprocess.run(cmd, stderr=None if debug else devnull, env=run_env)
    devnull.close()
    summary_file = mosdepth_prefix + ".mosdepth.summary.txt"
    quantized_bed = mosdepth_prefix + ".quantized.bed.gz"
    return summary_file, quantized_bed


def _read_quantized_bed(quantized_bed):
    """Read a mosdepth quantized BED file into a per-contig interval dict.

    Args:
        quantized_bed: Path to the .quantized.bed.gz file.

    Returns:
        Dict mapping contig name → list of (start, end, label) tuples
        in order of appearance.
    """
    data = {}
    opener = gzip.open if quantized_bed.endswith(".gz") else open
    with opener(quantized_bed, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            label = parts[3]
            if chrom not in data:
                data[chrom] = []
            data[chrom].append((start, end, label))
    return data


def _get_plot_prefix(input_file, report_file):
    """Derive the output path prefix for plot files.

    Strips .fasta/.fa/.fasta.gz/.fa.gz from the input basename and
    places the result in the same directory as the report file.

    Args:
        input_file: Path to the genome assembly FASTA.
        report_file: Path to the coverage report (used for output directory).

    Returns:
        String path prefix (no extension).
    """
    basename = os.path.basename(input_file)
    for ext in (".fasta.gz", ".fa.gz", ".fasta", ".fa"):
        if basename.endswith(ext):
            basename = basename[: -len(ext)]
            break
    else:
        basename = os.path.splitext(basename)[0]
    report_dir = os.path.dirname(os.path.abspath(report_file))
    return os.path.join(report_dir, basename)


def _cov_colours_for_labels(labels):
    """Return a colour dict for the given label list.

    Known labels get their canonical colour; others get sequential fallbacks.

    Args:
        labels: List of coverage-class label strings.

    Returns:
        Dict mapping label → hex colour string.
    """
    colours = {}
    fallback_idx = 0
    for lbl in labels:
        if lbl in _COV_COLOURS:
            colours[lbl] = _COV_COLOURS[lbl]
        else:
            colours[lbl] = _FALLBACK_COLOURS[fallback_idx % len(_FALLBACK_COLOURS)]
            fallback_idx += 1
    return colours


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------


def _save_figure(fig, path, plot_format):
    dpi = 150 if plot_format == "png" else None
    fig.savefig(path, format=plot_format, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_coverage_heatmap(contig_data, scaffold_rows, labels, colours, plot_prefix, plot_format):
    """Write per-scaffold coverage-class heatmap.

    Each scaffold is one horizontal row; regions are coloured by coverage
    class.  Scaffolds are paginated (_SCAFFOLDS_PER_PAGE per page).

    Args:
        contig_data: Dict from _read_quantized_bed.
        scaffold_rows: List of contig dicts from parse_mosdepth_summary
            (used for scaffold ordering).
        labels: Ordered list of coverage-class label strings.
        colours: Dict mapping label → hex colour.
        plot_prefix: Output path prefix (no extension).
        plot_format: "pdf", "svg", or "png".
    """
    scaffolds = [r["chrom"] for r in scaffold_rows if r["chrom"] in contig_data]
    if not scaffolds:
        return

    legend_patches = [Patch(facecolor=colours.get(lbl, "#888888"), label=lbl) for lbl in labels]
    pages = [scaffolds[i : i + _SCAFFOLDS_PER_PAGE] for i in range(0, len(scaffolds), _SCAFFOLDS_PER_PAGE)]

    def _make_page(page_scaffolds, page_num, total_pages):
        n = len(page_scaffolds)
        fig, ax = plt.subplots(figsize=(14, max(4, n * 0.35 + 2)))
        for row_idx, contig in enumerate(page_scaffolds):
            for s, e, lbl in contig_data.get(contig, []):
                ax.broken_barh(
                    [(s, e - s)],
                    (row_idx - 0.4, 0.8),
                    facecolors=colours.get(lbl, "#888888"),
                    linewidth=0,
                )
        ax.set_yticks(range(n))
        ax.set_yticklabels(page_scaffolds, fontsize=max(5, min(9, 200 // n)))
        ax.set_xlabel("Genomic position (bp)")
        title = "Coverage class heatmap"
        if total_pages > 1:
            title += f"  [page {page_num}/{total_pages}]"
        ax.set_title(title)
        ax.legend(handles=legend_patches, loc="upper right", fontsize=8)
        fig.tight_layout()
        return fig

    if plot_format == "pdf":
        path = plot_prefix + ".depth_heatmap.pdf"
        with PdfPages(path) as pdf:
            for i, page in enumerate(pages):
                fig = _make_page(page, i + 1, len(pages))
                pdf.savefig(fig)
                plt.close(fig)
        status(f"Coverage heatmap written to: {path}")
    else:
        for i, page in enumerate(pages):
            suffix = f"_p{i + 1:02d}" if len(pages) > 1 else ""
            path = f"{plot_prefix}.depth_heatmap{suffix}.{plot_format}"
            fig = _make_page(page, i + 1, len(pages))
            _save_figure(fig, path, plot_format)
            status(f"Coverage heatmap written to: {path}")


def _plot_coverage_barplot(contig_data, scaffold_rows, labels, colours, plot_prefix, plot_format):
    """Write per-scaffold stacked bar chart of coverage-class proportions.

    Args:
        contig_data: Dict from _read_quantized_bed.
        scaffold_rows: List of contig dicts from parse_mosdepth_summary.
        labels: Ordered list of coverage-class label strings.
        colours: Dict mapping label → hex colour.
        plot_prefix: Output path prefix (no extension).
        plot_format: "pdf", "svg", or "png".
    """
    scaffolds = [r["chrom"] for r in scaffold_rows if r["chrom"] in contig_data]
    if not scaffolds:
        return

    # Compute per-scaffold percentage of each coverage class
    pcts = {}
    for contig in scaffolds:
        intervals = contig_data.get(contig, [])
        total_bp = sum(e - s for s, e, _ in intervals)
        label_bp = {lbl: 0 for lbl in labels}
        for s, e, lbl in intervals:
            if lbl in label_bp:
                label_bp[lbl] += e - s
        pcts[contig] = {lbl: (label_bp[lbl] / total_bp * 100 if total_bp else 0.0) for lbl in labels}

    legend_patches = [Patch(facecolor=colours.get(lbl, "#888888"), label=lbl) for lbl in labels]
    pages = [scaffolds[i : i + _SCAFFOLDS_PER_PAGE] for i in range(0, len(scaffolds), _SCAFFOLDS_PER_PAGE)]

    def _make_page(page_scaffolds, page_num, total_pages):
        n = len(page_scaffolds)
        fig, ax = plt.subplots(figsize=(12, max(4, n * 0.35 + 2)))
        # draw bottom-to-top so first scaffold appears at the top
        display_order = list(reversed(page_scaffolds))
        for row_idx, contig in enumerate(display_order):
            left = 0.0
            for lbl in labels:
                pct = pcts[contig].get(lbl, 0.0)
                if pct > 0:
                    ax.barh(row_idx, pct, left=left, color=colours.get(lbl, "#888888"))
                left += pct
        ax.set_yticks(range(n))
        ax.set_yticklabels(display_order, fontsize=max(5, min(9, 200 // n)))
        ax.set_xlabel("Percentage of scaffold (bp)")
        ax.set_xlim(0, 100)
        title = "Coverage class proportions"
        if total_pages > 1:
            title += f"  [page {page_num}/{total_pages}]"
        ax.set_title(title)
        ax.legend(handles=legend_patches, loc="lower right", fontsize=8)
        fig.tight_layout()
        return fig

    if plot_format == "pdf":
        path = plot_prefix + ".depth_barplot.pdf"
        with PdfPages(path) as pdf:
            for i, page in enumerate(pages):
                fig = _make_page(page, i + 1, len(pages))
                pdf.savefig(fig)
                plt.close(fig)
        status(f"Coverage barplot written to: {path}")
    else:
        for i, page in enumerate(pages):
            suffix = f"_p{i + 1:02d}" if len(pages) > 1 else ""
            path = f"{plot_prefix}.depth_barplot{suffix}.{plot_format}"
            fig = _make_page(page, i + 1, len(pages))
            _save_figure(fig, path, plot_format)
            status(f"Coverage barplot written to: {path}")


def _plot_depth_histogram(contig_rows, mean_depth, plot_prefix, plot_format):
    """Write histogram and boxplot of per-scaffold mean coverage depths.

    Args:
        contig_rows: List of contig dicts from parse_mosdepth_summary.
        mean_depth: Assembly-wide mean depth (float).
        plot_prefix: Output path prefix (no extension).
        plot_format: "pdf", "svg", or "png".
    """
    depths = [c["mean"] for c in contig_rows if c.get("mean", 0) > 0]
    if not depths:
        return

    # Log-spaced bins for the histogram
    min_d = max(min(depths), 0.01)
    max_d = max(depths)
    log_min = math.log10(min_d)
    log_max = math.log10(max_d * 1.05)
    n_bins = min(60, len(depths))
    step = (log_max - log_min) / n_bins if log_max > log_min else 1.0
    bins = [10 ** (log_min + i * step) for i in range(n_bins + 1)]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram with log x-axis
    ax1 = axes[0]
    ax1.hist(depths, bins=bins, color="#2166ac", alpha=0.75, edgecolor="none")
    ax1.set_xscale("log")
    if mean_depth > 0:
        ax1.axvline(
            mean_depth,
            color="red",
            linestyle="--",
            linewidth=1,
            label=f"Mean: {mean_depth:.1f}x",
        )
        ax1.legend(fontsize=9)
    ax1.set_xlabel("Mean depth of coverage (log10 scale)")
    ax1.set_ylabel("Number of scaffolds")
    ax1.set_title("Per-scaffold coverage distribution")

    # Boxplot with log y-axis
    ax2 = axes[1]
    ax2.boxplot(
        depths,
        vert=True,
        patch_artist=True,
        boxprops=dict(facecolor="#2166ac", alpha=0.7),
        medianprops=dict(color="red", linewidth=1.5),
        flierprops=dict(marker="o", markersize=3, alpha=0.5),
    )
    ax2.set_yscale("log")
    ax2.set_ylabel("Mean depth of coverage (log10 scale)")
    ax2.set_xticks([1])
    ax2.set_xticklabels(["All scaffolds"])
    ax2.set_title("Coverage depth distribution")

    fig.suptitle("Scaffold coverage depth summary", fontsize=12, fontweight="bold")
    fig.tight_layout()

    path = f"{plot_prefix}.depth_histogram.{plot_format}"
    _save_figure(fig, path, plot_format)
    status(f"Depth histogram written to: {path}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def run(parser, args):
    """Execute depth-of-coverage analysis for a genome assembly.

    Maps Illumina and/or long reads to the assembly, computes per-contig
    and whole-assembly depth via mosdepth, and writes a coverage report.
    Contigs with mean depth > (assembly_mean + 3 * SD) are flagged as
    possible contaminants or organellar sequences.

    When plotting is not disabled (--no-plot), mosdepth is run in quantized
    mode and three coverage plots are produced alongside the text report.

    Args:
        parser: ArgumentParser instance (kept for CLI signature compatibility).
        args: Namespace of parsed CLI arguments.
    """
    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    if not args.left and not args.longreads:
        status("ERROR: provide at least --left (Illumina) or --longreads")
        sys.exit(1)

    if not checkfile(args.input):
        status(f"ERROR: assembly file not found or empty: {args.input}")
        sys.exit(1)

    genome = os.path.abspath(args.input)
    reads_left = os.path.abspath(args.left) if args.left else None
    reads_right = os.path.abspath(args.right) if args.right else None
    longreads = os.path.abspath(args.longreads) if args.longreads else None

    for label, fpath in [("--left", reads_left), ("--right", reads_right), ("--longreads", longreads)]:
        if fpath and not checkfile(fpath):
            status(f"ERROR: read file not found or empty ({label}): {fpath}")
            sys.exit(1)

    # ------------------------------------------------------------------
    # Check required tools
    # ------------------------------------------------------------------
    required = {"samtools", "mosdepth"}
    if reads_left:
        required.add(args.aligner)
    if longreads:
        required.add("minimap2")
    for tool in sorted(required):
        if not which(tool):
            status(f"ERROR: required tool '{tool}' not found in PATH. " f"Install via conda: conda install -c bioconda {tool}")
            sys.exit(1)

    no_plot = getattr(args, "no_plot", False)
    plot_format = getattr(args, "plot_format", "pdf")
    quantize_str = getattr(args, "quantize", _DEFAULT_QUANTIZE)
    quantize_labels_str = getattr(args, "quantize_labels", None)

    if not no_plot and not HAS_MATPLOTLIB:
        status("WARNING: matplotlib not available — coverage plots will be skipped.")
        status("         Install with: conda install -c conda-forge matplotlib")
        no_plot = True

    # ------------------------------------------------------------------
    # Working directory
    # ------------------------------------------------------------------
    custom_workdir = bool(args.workdir)
    if not args.workdir:
        args.workdir = f"aaftf-depth_{str(uuid.uuid4())[:8]}"
    workdir = os.path.abspath(args.workdir)
    if not os.path.exists(workdir):
        os.mkdir(workdir)

    report_file = args.out

    # ------------------------------------------------------------------
    # Count input reads
    # ------------------------------------------------------------------
    status("Counting input reads...")
    illumina_counts = {}
    if reads_left:
        status(f"  Counting reads in {os.path.basename(reads_left)}")
        illumina_counts["left"] = count_fastq_reads(reads_left)
    if reads_right:
        status(f"  Counting reads in {os.path.basename(reads_right)}")
        illumina_counts["right"] = count_fastq_reads(reads_right)
    lr_count = 0
    if longreads:
        status(f"  Counting reads in {os.path.basename(longreads)}")
        lr_count = count_fastq_reads(longreads)

    # ------------------------------------------------------------------
    # Map reads
    # ------------------------------------------------------------------
    status("Mapping reads to assembly...")
    bam_illumina, bam_longreads, bam_combined = map_reads(
        genome,
        reads_left,
        reads_right,
        longreads,
        workdir,
        args.cpus,
        args.illumina_preset,
        args.longread_preset,
        args.aligner,
        args.debug,
    )

    if not bam_combined or not os.path.exists(bam_combined):
        status("ERROR: mapping produced no BAM file")
        sys.exit(1)

    # ------------------------------------------------------------------
    # samtools flagstat
    # ------------------------------------------------------------------
    status("Running samtools flagstat...")
    flagstat_illumina = run_flagstat(bam_illumina) if bam_illumina else None
    flagstat_longreads = run_flagstat(bam_longreads) if bam_longreads else None

    # ------------------------------------------------------------------
    # mosdepth (quantized when plotting is enabled, standard otherwise)
    # ------------------------------------------------------------------
    status("Running mosdepth...")
    quantized_bed = None
    labels = None
    colours = None
    if not no_plot:
        labels, env_dict = _parse_quantize_bins(quantize_str, quantize_labels_str)
        colours = _cov_colours_for_labels(labels)
        summary_file, quantized_bed = run_mosdepth_quantized(
            bam_combined,
            workdir,
            args.cpus,
            quantize_str,
            env_dict,
            debug=args.debug,
        )
    else:
        summary_file = run_mosdepth(
            bam_combined,
            workdir,
            args.cpus,
            debug=args.debug,
        )

    if not os.path.exists(summary_file):
        status(f"ERROR: mosdepth summary not produced: {summary_file}")
        sys.exit(1)

    total_row, contig_rows = parse_mosdepth_summary(summary_file)

    dist_file = summary_file.replace(".mosdepth.summary.txt", ".mosdepth.global.dist.txt")
    coverage_breadth = _coverage_breadth_from_dist(dist_file)

    # ------------------------------------------------------------------
    # Statistics for outlier detection
    # ------------------------------------------------------------------
    min_contig_len = getattr(args, "min_contig_len", 500)
    analysis_rows = [c for c in contig_rows if c.get("length", 0) >= min_contig_len]
    if analysis_rows:
        depths = [c["mean"] for c in analysis_rows]
        mean_depth = total_row["mean"] if total_row else (sum(depths) / len(depths))
        # Population SD is intentional: the contigs ARE the full assembly
        # (not a sample). Sample SD would systematically push the threshold
        # above any single outlier, defeating outlier detection on small sets.
        sd_depth = math.sqrt(sum((d - mean_depth) ** 2 for d in depths) / len(depths))
        threshold = mean_depth + 3.0 * sd_depth
    else:
        mean_depth = total_row["mean"] if total_row else 0.0
        sd_depth = 0.0
        threshold = 0.0

    contig_rows_sorted = sorted(contig_rows, key=lambda x: x["mean"], reverse=True)
    n_outliers = sum(1 for c in contig_rows if c["mean"] > threshold)

    # ------------------------------------------------------------------
    # Write report
    # ------------------------------------------------------------------
    status(f"Writing coverage report to {report_file}")
    with open(report_file, "w") as fout:
        sep = "=" * 65
        fout.write(sep + "\n")
        fout.write("AAFTF Coverage Statistics Report\n")
        fout.write(sep + "\n")
        fout.write(f"Assembly: {genome}\n\n")

        # --- Section 1: Read Input Summary ---
        fout.write("=== 1. Read Input Summary ===\n")
        if reads_left:
            n = illumina_counts.get("left", -1)
            fout.write(f"  Illumina left reads:  {reads_left}\n")
            fout.write(f"    Read count:         {n:,}\n" if n >= 0 else "    Read count:         unknown\n")
        if reads_right:
            n = illumina_counts.get("right", -1)
            fout.write(f"  Illumina right reads: {reads_right}\n")
            fout.write(f"    Read count:         {n:,}\n" if n >= 0 else "    Read count:         unknown\n")
        if longreads:
            fout.write(f"  Long reads:           {longreads}\n")
            fout.write(f"    Read count:         {lr_count:,}\n" if lr_count >= 0 else "    Read count:         unknown\n")

        if flagstat_illumina:
            fout.write("\n  Illumina alignment (samtools flagstat):\n")
            for line in flagstat_illumina.splitlines():
                fout.write(f"    {line}\n")
        if flagstat_longreads:
            fout.write("\n  Long-read alignment (samtools flagstat):\n")
            for line in flagstat_longreads.splitlines():
                fout.write(f"    {line}\n")
        fout.write("\n")

        # --- Section 2: Whole-Assembly Coverage ---
        fout.write("=== 2. Whole-Assembly Coverage ===\n")
        if total_row:
            fout.write(f"  Mean depth:            {mean_depth:.2f}x\n")
            fout.write(f"  Total assembly length: {total_row['length']:,} bp\n")
            if coverage_breadth is not None:
                fout.write(f"  Bases covered (>=1x):  " f"{coverage_breadth * 100:.2f}%\n")
        fout.write("\n")

        # --- Section 3: Per-Contig Coverage ---
        fout.write("=== 3. Per-Contig Coverage (sorted by depth, descending) ===\n")
        fout.write(f"  Assembly mean: {mean_depth:.2f}x   SD: {sd_depth:.2f}x\n")
        fout.write(f"  Outlier threshold (mean + 3*SD): {threshold:.2f}x  " f"({n_outliers} contigs flagged)\n")
        fout.write("  Flagged contigs are possible contaminants or organellar sequences.\n\n")

        cw = (40, 14, 12)
        header = f"  {'Contig':<{cw[0]}} " f"{'Length (bp)':>{cw[1]}} " f"{'Mean Depth':>{cw[2]}}  Flag\n"
        fout.write(header)
        fout.write("  " + "-" * (sum(cw) + 10) + "\n")
        for c in contig_rows_sorted:
            flag = "** OUTLIER (possible contaminant/organelle)" if c["mean"] > threshold else ""
            fout.write(f"  {c['chrom']:<{cw[0]}} " f"{c['length']:>{cw[1]},} " f"{c['mean']:>{cw[2]}.2f}  {flag}\n")

    status(f"Coverage report written to: {report_file}")

    # ------------------------------------------------------------------
    # Coverage plots
    # ------------------------------------------------------------------
    if quantized_bed and os.path.exists(quantized_bed):
        status("Reading quantized coverage BED...")
        contig_data = _read_quantized_bed(quantized_bed)
        plot_prefix = _get_plot_prefix(args.input, args.out)
        status("Generating coverage plots...")
        _plot_coverage_heatmap(contig_data, contig_rows_sorted, labels, colours, plot_prefix, plot_format)
        _plot_coverage_barplot(contig_data, contig_rows_sorted, labels, colours, plot_prefix, plot_format)
        _plot_depth_histogram(contig_rows, mean_depth, plot_prefix, plot_format)

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------
    if not args.debug and not custom_workdir:
        SafeRemove(workdir)

    if not args.pipe:
        status(f"Your next command might be:\n" f"\tAAFTF assess -i {args.input}\n")
