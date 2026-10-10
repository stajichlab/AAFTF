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

import logging
import math
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

from aaftf.utility import align_to_sorted_bam, check_file, cleanup_workdir, close_step_log, count_fastq, make_workdir, open_maybe_gz, print_cmd, require_tools, run_cmd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Patch

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

__all__ = ["HAS_MATPLOTLIB", "run", "map_reads", "run_flagstat", "run_mosdepth", "parse_mosdepth_summary"]


logger = logging.getLogger(__name__)

_DEFAULT_QUANTIZE = "0:1:4:100:200:"

_COV_COLOURS = {
    "NO_COVERAGE": "#000080",
    "LOW_COVERAGE": "#006fff",
    "CALLABLE": "#54ffaa",
    "HIGH_COVERAGE": "#ff6c00",
    "VERY_HIGH_COVERAGE": "#800000",
}

_DEFAULT_LABELS = list(_COV_COLOURS)

_SCAFFOLDS_PER_PAGE = 50

_HEATMAP_MAX_LENGTH_RATIO = 10  # max longest/shortest contig length ratio per page


def run(
    input: str,
    out: str = "coverage_stats.txt",
    read1: str | None = None,
    read2: str | None = None,
    longreads: str | None = None,
    longread_preset: str | None = None,
    aligner: str = "minimap2",
    cpus: int = 1,
    workdir: str | None = None,
    debug: bool = False,
    pipe: bool = False,
    min_contig_len: int = 500,
    no_plot: bool = False,
    plot_format: str = "pdf",
    **kwargs: Any,
) -> None:
    """Execute depth-of-coverage analysis for a genome assembly.

    Maps Illumina and/or long reads to the assembly, computes per-contig
    and whole-assembly depth via mosdepth, and writes a coverage report.
    Contigs with mean depth > (assembly_mean + 3 * SD) are flagged as
    possible contaminants or organellar sequences.

    When plotting is not disabled (--no-plot), mosdepth is run in quantized
    mode and three coverage plots are produced alongside the text report.

    Args:
        input: Path to the genome assembly FASTA.
        out: Path of the coverage report to write.
        read1: Forward (or single-end) Illumina FASTQ, or None.
        read2: Reverse Illumina FASTQ, or None.
        longreads: Long-read FASTQ, or None.
        longread_preset: minimap2 preset for long reads (``map-ont``, ``map-pb`` or ``map-hifi``);
            required when ``longreads`` is given.
        aligner: Aligner for Illumina reads (``minimap2`` or ``bwa``).
        cpus: Number of CPU threads.
        workdir: Working directory for intermediate files; a default is created when None.
        debug: If True, show subprocess stderr and keep the working directory.
        min_contig_len: Minimum contig length included in the outlier statistics.
        no_plot: If True, skip quantized mosdepth and plot generation.
        plot_format: Plot file format (``pdf``, ``svg`` or ``png``).
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``); ignored.

    Raises:
        ValueError: If neither ``read1`` nor ``longreads`` is given, or ``longreads`` lacks a preset.
        FileNotFoundError: If the assembly or a read file is missing or empty.
        RuntimeError: If mapping produces no BAM or mosdepth produces no summary file.
    """
    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    if not read1 and not longreads:
        raise ValueError("provide at least --read1 (Illumina) or --longreads")

    if longreads and not longread_preset:
        raise ValueError("--longread_preset is required when --longreads is provided (map-ont, map-pb, or map-hifi)")

    if not check_file(input):
        raise FileNotFoundError(f"assembly file not found or empty: {input}")

    genome = str(Path(input).resolve())
    read1 = str(Path(read1).resolve()) if read1 else None
    read2 = str(Path(read2).resolve()) if read2 else None
    longreads = str(Path(longreads).resolve()) if longreads else None

    for label, fpath in [("--read1", read1), ("--read2", read2), ("--longreads", longreads)]:
        if fpath and not check_file(fpath):
            raise FileNotFoundError(f"read file not found or empty ({label}): {fpath}")

    # ------------------------------------------------------------------
    # Check required tools
    # ------------------------------------------------------------------
    required = {"samtools", "mosdepth"}
    if read1:
        required.add(aligner)
    if longreads:
        required.add("minimap2")
    require_tools(sorted(required))

    if not no_plot and not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available — coverage plots will be skipped.")
        logger.info("Install with: conda install -c conda-forge matplotlib")
        no_plot = True

    # ------------------------------------------------------------------
    # Working directory
    # ------------------------------------------------------------------
    workdir, custom_workdir = make_workdir(workdir, "depth")
    try:
        workdir = str(Path(workdir).resolve())

        report_file = out

        # ------------------------------------------------------------------
        # Count input reads
        # ------------------------------------------------------------------
        logger.info("Counting input reads...")
        read_counts = {}  # -1 marks a file that could not be read; reported as "unknown"
        for key, fastq in (("read1", read1), ("read2", read2), ("long", longreads)):
            if fastq:
                logger.info(f"Counting reads in {Path(fastq).name}")
                try:
                    read_counts[key] = count_fastq(fastq)
                except (OSError, subprocess.CalledProcessError):
                    read_counts[key] = -1

        # ------------------------------------------------------------------
        # Map reads
        # ------------------------------------------------------------------
        logger.info("Mapping reads to assembly...")
        bam_illumina, bam_longreads, bam_combined = map_reads(
            genome,
            read1,
            read2,
            longreads,
            workdir,
            cpus,
            "sr",
            longread_preset,
            aligner,
            debug,
        )

        if not bam_combined or not Path(bam_combined).exists():
            raise RuntimeError("mapping produced no BAM file")

        # ------------------------------------------------------------------
        # samtools flagstat
        # ------------------------------------------------------------------
        logger.info("Running samtools flagstat...")
        flagstat_illumina = run_flagstat(bam_illumina) if bam_illumina else None
        flagstat_longreads = run_flagstat(bam_longreads) if bam_longreads else None

        # ------------------------------------------------------------------
        # mosdepth (quantized when plotting is enabled, standard otherwise)
        # ------------------------------------------------------------------
        logger.info("Running mosdepth...")
        labels = None
        colors = None
        if not no_plot:
            labels, colors = _parse_quantize_bins()
        summary_file, quantized_bed = run_mosdepth(
            bam_combined,
            workdir,
            cpus,
            labels=labels,
            debug=debug,
        )

        if not Path(summary_file).exists():
            raise RuntimeError(f"mosdepth summary not produced: {summary_file}")

        total_row, contig_rows = parse_mosdepth_summary(summary_file)

        coverage_breadth = _coverage_breadth_from_dist(summary_file.with_name(summary_file.name.removesuffix(".mosdepth.summary.txt") + ".mosdepth.global.dist.txt"))

        # ------------------------------------------------------------------
        # Statistics for outlier detection
        # ------------------------------------------------------------------
        stats = _depth_outliers(contig_rows, total_row, min_contig_len)
        mosdepth_mean_depth = stats["mosdepth_mean_depth"]
        contig_arith_mean = stats["contig_arith_mean"]
        mean_depth = stats["mean_depth"]
        sd_depth = stats["sd_depth"]
        threshold_2sd = stats["threshold_2sd"]
        threshold_3sd = stats["threshold_3sd"]
        n_outliers_2sd = stats["n_outliers_2sd"]
        n_outliers_3sd = stats["n_outliers_3sd"]
        contig_depth_sorted = sorted(contig_rows, key=lambda x: x["mean"], reverse=True)
        contig_length_sorted = sorted(contig_rows, key=lambda x: x["length"], reverse=True)

        # ------------------------------------------------------------------
        # Write report
        # ------------------------------------------------------------------
        logger.info(f"Writing coverage report to {report_file}")
        with open(report_file, "w") as fout:
            sep = "=" * 65
            fout.write(sep + "\n")
            fout.write("AAFTF Coverage Statistics Report\n")
            fout.write(sep + "\n")
            fout.write(f"Assembly: {genome}\n\n")

            # --- Section 1: Read Input Summary ---
            fout.write("=== 1. Read Input Summary ===\n")
            if read1:
                n = read_counts["read1"]
                fout.write(f"  Illumina read 1:  {read1}\n")
                fout.write(f"    Read count:         {n:,}\n" if n >= 0 else "    Read count:         unknown\n")
            if read2:
                n = read_counts["read2"]
                fout.write(f"  Illumina read 2:  {read2}\n")
                fout.write(f"    Read count:         {n:,}\n" if n >= 0 else "    Read count:         unknown\n")
            if longreads:
                fout.write(f"  Long reads:           {longreads}\n")
                n = read_counts["long"]
                fout.write(f"    Read count:         {n:,}\n" if n >= 0 else "    Read count:         unknown\n")

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
                fout.write(f"  Mean depth (mosdepth global, length-weighted): {mosdepth_mean_depth:.2f}x\n")
                fout.write(f"  Mean depth (per-contig arithmetic mean):        {contig_arith_mean:.2f}x\n")
                fout.write(f"  Total assembly length: {total_row['length']:,} bp\n")
                if coverage_breadth is not None:
                    fout.write(f"  Bases covered (>=1x):  {coverage_breadth * 100:.2f}%\n")
            else:
                fout.write(f"  Mean depth (per-contig arithmetic mean): {contig_arith_mean:.2f}x\n")
            fout.write("\n")

            # --- Section 3: Per-Contig Coverage ---
            fout.write("=== 3. Per-Contig Coverage (sorted by depth, descending) ===\n")
            fout.write(f"  Assembly mean: {mean_depth:.2f}x   SD: {sd_depth:.2f}x\n")
            fout.write(f"  Elevated threshold  (mean + 2*SD): {threshold_2sd:.2f}x  ({n_outliers_2sd} contigs)\n")
            fout.write(f"  Outlier threshold   (mean + 3*SD): {threshold_3sd:.2f}x  ({n_outliers_3sd} contigs)\n")
            fout.write("  OUTLIER contigs are likely contaminants or organellar sequences.\n")
            fout.write("  ELEVATED contigs are candidates worth inspecting (2–3 SD above mean).\n\n")

            cw = (40, 14, 12)
            header = f"  {'Contig':<{cw[0]}} {'Length (bp)':>{cw[1]}} {'Mean Depth':>{cw[2]}}  Flag\n"
            fout.write(header)
            fout.write("  " + "-" * (sum(cw) + 10) + "\n")
            for c in contig_depth_sorted:
                depth = c["mean"]
                if depth > threshold_3sd:
                    flag = "** OUTLIER (possible contaminant/organelle)"
                elif depth > threshold_2sd:
                    flag = "   ELEVATED (inspect — 2–3 SD above mean)"
                else:
                    flag = ""
                fout.write(f"  {c['chrom']:<{cw[0]}} {c['length']:>{cw[1]},} {c['mean']:>{cw[2]}.2f}  {flag}\n")

        logger.info(f"Coverage report written to: {report_file}")

        # ------------------------------------------------------------------
        # Coverage plots
        # ------------------------------------------------------------------
        if labels and colors and quantized_bed and Path(quantized_bed).exists():
            logger.info("Reading quantized coverage BED...")
            contig_data = _read_quantized_bed(quantized_bed)
            plot_prefix = _get_plot_prefix(input, out)
            logger.info("Generating coverage plots...")
            _plot_coverage_heatmap(contig_data, contig_length_sorted, labels, colors, plot_prefix, plot_format)
            _plot_coverage_barplot(contig_data, contig_depth_sorted, labels, colors, plot_prefix, plot_format)
            _plot_depth_histogram(contig_rows, mean_depth, plot_prefix, plot_format)

        # ------------------------------------------------------------------
        # Cleanup
        # ------------------------------------------------------------------
        cleanup_workdir(workdir, debug, custom_workdir)

        if not pipe:
            logger.info(f"Your next command might be:\nAAFTF assess -i {input}")
    finally:
        close_step_log()


def map_reads(
    genome: str,
    read1: str | None,
    read2: str | None,
    longreads: str | None,
    workdir: str,
    cpus: int,
    illumina_preset: str,
    longread_preset: str | None,
    aligner: str,
    debug: bool,
) -> tuple[str | None, str | None, str | None]:
    """Map reads to the genome assembly and produce sorted, indexed BAM files.

    Illumina reads are mapped with minimap2 (default) or bwa mem.
    Long reads are always mapped with minimap2.
    When both read types are provided, the two BAM files are merged.

    Args:
        genome: Absolute path to genome assembly FASTA.
        read1: Path to forward Illumina FASTQ (or None).
        read2: Path to reverse Illumina FASTQ (or None).
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
        mosdepth (merged when both types are present); None if no reads were given.

    Raises:
        RuntimeError: If ``bwa index``, or ``samtools merge``/``index`` of the combined BAM, fails.
        ValueError: If ``longreads`` is given without ``longread_preset``.
    """
    bam_illumina = None
    bam_longreads = None

    # --- Illumina reads ---
    if read1:
        bam_illumina = str(Path(workdir, "illumina.sorted.bam"))
        if aligner == "bwa":
            # Copy genome to workdir so BWA index files stay contained
            genome_local = str(Path(workdir, Path(genome).name))
            shutil.copyfile(genome, genome_local)
            bwa_index_cmd = ["bwa", "index", genome_local]
            ret = run_cmd(bwa_index_cmd, debug)
            if ret.returncode != 0:
                raise RuntimeError("bwa index failed")
            read_group = r"@RG\tID:illumina\tSM:illumina\tPL:illumina"
            map_cmd = ["bwa", "mem", "-t", str(cpus), "-R", read_group, genome_local, read1]
            if read2:
                map_cmd.append(read2)
        else:
            map_cmd = ["minimap2", "-ax", illumina_preset, "-t", str(cpus), genome, read1]
            if read2:
                map_cmd.append(read2)

        align_to_sorted_bam(map_cmd, bam_illumina, cpus, debug=debug)

    # --- Long reads ---
    if longreads:
        if not longread_preset:
            raise ValueError("--longread_preset is required when --longreads is provided (map-ont, map-pb, or map-hifi)")
        bam_longreads = str(Path(workdir, "longreads.sorted.bam"))
        map_cmd = ["minimap2", "-ax", longread_preset, "-t", str(cpus), genome, longreads]
        align_to_sorted_bam(map_cmd, bam_longreads, cpus, debug=debug)

    # --- Combine ---
    bam_combined: str | None
    if bam_illumina and bam_longreads:
        bam_combined = str(Path(workdir, "combined.sorted.bam"))
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
        if run_cmd(merge_cmd, debug).returncode != 0 or run_cmd(["samtools", "index", bam_combined], debug).returncode != 0:
            raise RuntimeError(f"samtools merge/index failed for {bam_combined}")
    elif bam_illumina:
        bam_combined = bam_illumina
    else:
        bam_combined = bam_longreads

    return bam_illumina, bam_longreads, bam_combined


def run_flagstat(bam_file: str) -> str:
    """Run samtools flagstat on a BAM file.

    Args:
        bam_file: Path to indexed BAM file.

    Returns:
        String containing the flagstat output lines.

    Raises:
        RuntimeError: If samtools flagstat fails.
    """
    cmd = ["samtools", "flagstat", bam_file]
    print_cmd(cmd)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"samtools flagstat failed for {bam_file}: {result.stderr.strip()}")
    return result.stdout


def run_mosdepth(
    bam_file: str,
    workdir: str,
    cpus: int,
    labels: list[str] | None = None,
    quantize_str: str = _DEFAULT_QUANTIZE,
    prefix: str = "coverage",
    debug: bool = False,
) -> tuple[Path, Path | None]:
    """Run mosdepth to calculate per-contig depth statistics.

    When `labels` is given, mosdepth is run in quantized mode: a quantized
    BED file is also produced, showing which coverage class each genomic
    interval belongs to (labels are set as MOSDEPTH_Q? env vars for the BED).

    Args:
        bam_file: Absolute path to the sorted, indexed BAM file.
        workdir: Directory for mosdepth output files.
        cpus: Number of CPU threads.
        labels: Optional list of bin label strings; enables quantized mode.
        quantize_str: Colon-separated quantize boundaries, e.g. "0:1:4:100:200:".
            Only used when `labels` is given.
        prefix: Filename prefix for mosdepth outputs.
        debug: If True, show mosdepth stderr.

    Returns:
        Tuple (summary_file, quantized_bed) — Path objects for the mosdepth
        summary text file and the quantized BED file. quantized_bed is None
        unless `labels` was given.
    """
    mosdepth_prefix = Path(workdir).resolve() / prefix
    cmd = ["mosdepth", "-n", "--threads", str(cpus)]
    run_env = None
    quantized_bed = None
    if labels is not None:
        run_env = os.environ.copy()
        run_env.update({f"MOSDEPTH_Q{i}": lbl for i, lbl in enumerate(labels)})
        cmd += ["--quantize", quantize_str]
        quantized_bed = mosdepth_prefix.with_name(mosdepth_prefix.name + ".quantized.bed.gz")
    cmd += [str(mosdepth_prefix), str(Path(bam_file).resolve())]
    run_cmd(cmd, debug, env=run_env)
    summary_file = mosdepth_prefix.with_name(mosdepth_prefix.name + ".mosdepth.summary.txt")
    return summary_file, quantized_bed


def parse_mosdepth_summary(
    summary_file: str | Path,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
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


def _depth_outliers(contig_rows: list[dict[str, Any]], total_row: dict[str, Any] | None, min_contig_len: int) -> dict[str, Any]:
    """Compute assembly depth statistics and outlier counts.

    Args:
        contig_rows: Per-contig mosdepth rows with ``length`` and ``mean`` keys.
        total_row: The mosdepth ``total`` row, or None when absent.
        min_contig_len: Contigs shorter than this are excluded from the
            mean/SD statistics (but still counted as outliers).

    Returns:
        Dict with ``mosdepth_mean_depth`` (float or None), ``contig_arith_mean``,
        ``mean_depth``, ``sd_depth`` (population SD), ``threshold_2sd``,
        ``threshold_3sd``, ``n_outliers_2sd`` (contigs in (2SD, 3SD]) and
        ``n_outliers_3sd`` (contigs above 3SD).
    """
    analysis_rows = [c for c in contig_rows if c.get("length", 0) >= min_contig_len]
    # mosdepth global mean: bases covered / assembly length (length-weighted)
    mosdepth_mean_depth = total_row["mean"] if total_row else None
    if analysis_rows:
        depths = [c["mean"] for c in analysis_rows]
        # Per-contig arithmetic mean (unweighted)
        contig_arith_mean = sum(depths) / len(depths)
        # Use mosdepth global mean for outlier threshold when available; it is
        # the more accurate estimate because it weights by contig length.
        mean_depth = mosdepth_mean_depth if mosdepth_mean_depth is not None else contig_arith_mean
        # Population SD is intentional: the contigs ARE the full assembly
        # (not a sample). Sample SD would systematically push the threshold
        # above any single outlier, defeating outlier detection on small sets.
        sd_depth = math.sqrt(sum((d - mean_depth) ** 2 for d in depths) / len(depths))
        threshold_2sd = mean_depth + 2.0 * sd_depth
        threshold_3sd = mean_depth + 3.0 * sd_depth
    else:
        contig_arith_mean = mosdepth_mean_depth if mosdepth_mean_depth is not None else 0.0
        mean_depth = contig_arith_mean
        sd_depth = 0.0
        threshold_2sd = 0.0
        threshold_3sd = 0.0

    n_outliers_2sd = sum(1 for c in contig_rows if threshold_2sd < c["mean"] <= threshold_3sd)
    n_outliers_3sd = sum(1 for c in contig_rows if c["mean"] > threshold_3sd)

    return {
        "mosdepth_mean_depth": mosdepth_mean_depth,
        "contig_arith_mean": contig_arith_mean,
        "mean_depth": mean_depth,
        "sd_depth": sd_depth,
        "threshold_2sd": threshold_2sd,
        "threshold_3sd": threshold_3sd,
        "n_outliers_2sd": n_outliers_2sd,
        "n_outliers_3sd": n_outliers_3sd,
    }


def _parse_quantize_bins() -> tuple[list[str], dict[str, str]]:
    """Return the default mosdepth quantize bin labels and their colors.

    Known labels get their canonical color; others fall back to grey.

    Returns:
        Tuple (labels, colors) where labels is a list of strings and colors
        is a dict mapping label -> hex color string.
    """
    labels = list(_DEFAULT_LABELS)
    colors = {lbl: _COV_COLOURS.get(lbl, "#888888") for lbl in labels}
    return labels, colors


def _coverage_breadth_from_dist(dist_file: str | Path) -> float | None:
    """Read coverage breadth at >= 1x from a mosdepth global distribution file.

    Any error while reading the file is swallowed and reported as None.

    Args:
        dist_file: Path to *.mosdepth.global.dist.txt.

    Returns:
        Float fraction (0–1) of bases with depth >= 1, or None if unavailable.
    """
    if not Path(dist_file).exists():
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


def _read_quantized_bed(quantized_bed: str | Path) -> dict[str, list[tuple[int, int, str]]]:
    """Read a mosdepth quantized BED file into a per-contig interval dict.

    Args:
        quantized_bed: Path to the .quantized.bed.gz file.

    Returns:
        Dict mapping contig name → list of (start, end, label) tuples
        in order of appearance.
    """
    data: dict[str, list[tuple[int, int, str]]] = {}
    with open_maybe_gz(quantized_bed) as fh:
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


def _get_plot_prefix(input_file: str | Path, report_file: str | Path) -> str:
    """Derive the output path prefix for plot files.

    Strips .fasta/.fa/.fasta.gz/.fa.gz from the input basename (otherwise the
    last suffix) and places the result in the same directory as the report file.

    Args:
        input_file: Path to the genome assembly FASTA.
        report_file: Path to the coverage report (used for output directory).

    Returns:
        String path prefix (no extension).
    """
    basename = Path(input_file).name
    for ext in (".fasta.gz", ".fa.gz", ".fasta", ".fa"):
        if basename.endswith(ext):
            basename = basename[: -len(ext)]
            break
    else:
        basename = Path(basename).stem
    report_dir = str(Path(report_file).resolve().parent)
    return str(Path(report_dir, basename))


def _plot_coverage_heatmap(
    contig_data: dict[str, list[tuple[int, int, str]]],
    scaffold_rows: list[dict[str, Any]],
    labels: list[str],
    colors: dict[str, str],
    plot_prefix: str,
    plot_format: str,
) -> None:
    """Write per-scaffold coverage-class heatmap.

    Each scaffold is one horizontal row; regions are colored by coverage class. Scaffolds are paginated by length ratio (see
    _paginate_by_length_ratio) so contigs of wildly different sizes don't share a page and squash each other on the shared bp
    x-axis.

    Args:
        contig_data: Dict from _read_quantized_bed.
        scaffold_rows: List of contig dicts from parse_mosdepth_summary,
            sorted by length descending.
        labels: Ordered list of coverage-class label strings.
        colors: Dict mapping label → hex color.
        plot_prefix: Output path prefix (no extension).
        plot_format: "pdf", "svg", or "png".
    """
    present_rows = [r for r in scaffold_rows if r["chrom"] in contig_data]
    if not present_rows:
        return

    legend_patches = [Patch(facecolor=colors.get(lbl, "#888888"), label=lbl) for lbl in labels]
    pages = _paginate_by_length_ratio(present_rows)
    # Fixed across pages (derived from the largest page) so labels are the
    # same size everywhere, regardless of how many rows a given page has.
    ytick_fontsize = max(5, min(9, 200 // max(len(p) for p in pages)))

    def _make_page(scaffold_idx_list: list[str], page_idx: int, total_pages: int) -> Any:
        """Build one heatmap page figure.

        Args:
            scaffold_idx_list: Contig names to draw on this page, one per row.
            page_idx: 1-based page number (shown in the title).
            total_pages: Total number of pages.

        Returns:
            The matplotlib Figure for this page.
        """
        n = len(scaffold_idx_list)
        fig_h = max(4, n * 0.35 + 2)
        fig, ax = plt.subplots(figsize=(14, fig_h))
        for row_idx, contig in enumerate(scaffold_idx_list):
            for s, e, lbl in contig_data.get(contig, []):
                ax.broken_barh(
                    [(s, e - s)],
                    (row_idx - 0.4, 0.8),
                    facecolors=colors.get(lbl, "#888888"),
                    linewidth=0,
                )
        ax.set_ylim(-0.5, n - 0.5)
        ax.set_yticks(range(n))
        ax.set_yticklabels(scaffold_idx_list, fontsize=ytick_fontsize)
        ax.invert_yaxis()
        ax.set_xlabel("Genomic position (bp)")
        title = "Coverage class heatmap"
        if total_pages > 1:
            title += f"  [page {page_idx}/{total_pages}]"
        # Fixed inch offsets from the top of the figure keep the title/legend
        # at a consistent visual position regardless of per-page figure height.
        title_in, legend_in, axes_top_in = 0.3, 0.65, 1.0
        fig.suptitle(title, y=1 - title_in / fig_h)
        fig.legend(
            handles=legend_patches,
            loc="upper center",
            bbox_to_anchor=(0.5, 1 - legend_in / fig_h),
            bbox_transform=fig.transFigure,
            ncol=5,
            fontsize=8,
        )
        fig.tight_layout(rect=(0, 0, 1, 1 - axes_top_in / fig_h))
        return fig

    if plot_format == "pdf":
        path = plot_prefix + ".depth_heatmap.pdf"
        with PdfPages(path) as pdf:
            for page_idx, scaffold_idx_list in enumerate(pages):
                fig = _make_page(scaffold_idx_list, page_idx + 1, len(pages))
                pdf.savefig(fig)
                plt.close(fig)
        logger.info(f"Coverage heatmap written to: {path}")
    else:
        for i, page in enumerate(pages):
            suffix = f"_p{i + 1:02d}" if len(pages) > 1 else ""
            path = f"{plot_prefix}.depth_heatmap{suffix}.{plot_format}"
            fig = _make_page(page, i + 1, len(pages))
            _save_figure(fig, path, plot_format)
            logger.info(f"Coverage heatmap written to: {path}")


def _paginate_by_length_ratio(
    scaffold_rows: list[dict[str, Any]],
    max_ratio: float = _HEATMAP_MAX_LENGTH_RATIO,
    max_per_page: int = _SCAFFOLDS_PER_PAGE,
) -> list[list[str]]:
    """Group scaffolds (pre-sorted by length, descending) into pages.

    Starts a new page once the running page's longest/shortest length
    ratio would exceed max_ratio, or once max_per_page rows is reached.

    Args:
        scaffold_rows: List of contig dicts (with 'chrom' and 'length'),
            sorted by length descending.
        max_ratio: Max allowed ratio between the longest and shortest
            contig length within one page.
        max_per_page: Hard cap on rows per page.

    Returns:
        List of pages, each a list of chrom name strings.
    """
    pages: list[list[str]] = []
    current: list[str] = []
    page_max_len = 0
    for row in scaffold_rows:
        length = row.get("length") or 1
        if current and (length < page_max_len / max_ratio or len(current) >= max_per_page):
            pages.append(current)
            current = []
        if not current:
            page_max_len = length
        current.append(row["chrom"])
    if current:
        pages.append(current)
    return pages


def _save_figure(fig: Any, path: str, plot_format: str) -> None:
    """Save a matplotlib figure to disk and close it.

    PNG output is rendered at 150 dpi; other formats use the matplotlib default.

    Args:
        fig: matplotlib Figure to save.
        path: Output file path.
        plot_format: File format passed to ``savefig`` (``pdf``, ``svg`` or ``png``).
    """
    dpi = 150 if plot_format == "png" else None
    fig.savefig(path, format=plot_format, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_coverage_barplot(
    contig_data: dict[str, list[tuple[int, int, str]]],
    scaffold_rows: list[dict[str, Any]],
    labels: list[str],
    colors: dict[str, str],
    plot_prefix: str,
    plot_format: str,
) -> None:
    """Write per-scaffold stacked bar chart of coverage-class proportions.

    Args:
        contig_data: Dict from _read_quantized_bed.
        scaffold_rows: List of contig dicts from parse_mosdepth_summary.
        labels: Ordered list of coverage-class label strings.
        colors: Dict mapping label → hex color.
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

    legend_patches = [Patch(facecolor=colors.get(lbl, "#888888"), label=lbl) for lbl in labels]
    pages = [scaffolds[i : i + _SCAFFOLDS_PER_PAGE] for i in range(0, len(scaffolds), _SCAFFOLDS_PER_PAGE)]

    def _make_page(page_scaffolds: list[str], page_num: int, total_pages: int) -> Any:
        """Build one stacked-bar page figure.

        Args:
            page_scaffolds: Contig names to draw on this page, one per bar.
            page_num: 1-based page number (shown in the title).
            total_pages: Total number of pages.

        Returns:
            The matplotlib Figure for this page.
        """
        n = len(page_scaffolds)
        fig_h = max(4, n * 0.35 + 2)
        fig, ax = plt.subplots(figsize=(12, fig_h))
        for row_idx, contig in enumerate(page_scaffolds):
            left = 0.0
            for lbl in labels:
                pct = pcts[contig].get(lbl, 0.0)
                if pct > 0:
                    ax.barh(row_idx, pct, left=left, color=colors.get(lbl, "#888888"))
                left += pct
        ax.set_ylim(-0.5, n - 0.5)
        ax.set_yticks(range(n))
        ax.set_yticklabels(page_scaffolds, fontsize=max(5, min(9, 200 // _SCAFFOLDS_PER_PAGE)))
        ax.invert_yaxis()
        ax.set_xlabel("Percentage of scaffold (bp)")
        ax.set_xlim(0, 100)
        title = "Coverage class proportions"
        if total_pages > 1:
            title += f"  [page {page_num}/{total_pages}]"
        # Fixed inch offsets from the top of the figure keep the title/legend
        # at a consistent visual position regardless of per-page figure height.
        title_in, legend_in, axes_top_in = 0.3, 0.65, 1.0
        fig.suptitle(title, y=1 - title_in / fig_h)
        fig.legend(
            handles=legend_patches,
            loc="upper center",
            bbox_to_anchor=(0.5, 1 - legend_in / fig_h),
            bbox_transform=fig.transFigure,
            ncol=5,
            fontsize=8,
        )
        fig.tight_layout(rect=(0, 0, 1, 1 - axes_top_in / fig_h))
        return fig

    if plot_format == "pdf":
        path = plot_prefix + ".depth_barplot.pdf"
        with PdfPages(path) as pdf:
            for i, page in enumerate(pages):
                fig = _make_page(page, i + 1, len(pages))
                pdf.savefig(fig)
                plt.close(fig)
        logger.info(f"Coverage barplot written to: {path}")
    else:
        for i, page in enumerate(pages):
            suffix = f"_p{i + 1:02d}" if len(pages) > 1 else ""
            path = f"{plot_prefix}.depth_barplot{suffix}.{plot_format}"
            fig = _make_page(page, i + 1, len(pages))
            _save_figure(fig, path, plot_format)
            logger.info(f"Coverage barplot written to: {path}")


def _plot_depth_histogram(contig_rows: list[dict[str, Any]], mean_depth: float, plot_prefix: str, plot_format: str) -> None:
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
        depths,  # vertical is the default orientation (vert= is deprecated in matplotlib >= 3.10)
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
    logger.info(f"Depth histogram written to: {path}")
