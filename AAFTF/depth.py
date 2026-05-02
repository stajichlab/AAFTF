"""Calculate depth of coverage of reads mapped to a genome assembly.

Maps Illumina paired-end and/or long reads to the genome assembly using
minimap2 (or bwa for Illumina reads), then runs mosdepth to compute
per-contig and whole-assembly coverage statistics.  Contigs with mean
depth more than 3 standard deviations above the assembly mean are
flagged as possible contaminants or organelles.
"""

import gzip
import math
import os
import shutil
import subprocess
import sys
import uuid

from AAFTF.utility import SafeRemove, checkfile, printCMD, status, which

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
# Main entry point
# ---------------------------------------------------------------------------


def run(parser, args):
    """Execute depth-of-coverage analysis for a genome assembly.

    Maps Illumina and/or long reads to the assembly, computes per-contig
    and whole-assembly depth via mosdepth, and writes a coverage report.
    Contigs with mean depth > (assembly_mean + 3 * SD) are flagged as
    possible contaminants or organellar sequences.

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
    # mosdepth
    # ------------------------------------------------------------------
    status("Running mosdepth...")
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
    # Cleanup
    # ------------------------------------------------------------------
    if not args.debug and not custom_workdir:
        SafeRemove(workdir)

    if not args.pipe:
        status(f"Your next command might be:\n" f"\tAAFTF assess -i {args.input}\n")
