"""Remove contigs with low read coverage relative to the bulk of the assembly."""

import logging
import shutil
from pathlib import Path
from typing import Any

from Bio import SeqIO

from aaftf.utility import align_to_sorted_bam, calc_nx, cleanup_workdir, close_step_log, execute, fasta_stats, filter_fasta, make_workdir, next_step_name, run_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(
    input: str,
    outfile: str,
    read1: str,
    read2: str | None = None,
    workdir: str | None = None,
    cpus: int = 1,
    mincovpct: int = 5,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Map reads to the assembly with BWA and drop contigs with low coverage.

    The baseline is the average coverage of the contigs at least N50 long; contigs whose mean
    coverage is at most ``mincovpct`` percent of it are dropped.

    Args:
        input: Assembly FASTA to screen.
        outfile: Output FASTA of retained contigs.
        read1: Forward (or single-end) reads FASTQ.
        read2: Reverse reads FASTQ for paired-end data.
        workdir: Working directory; a temporary one is created when None.
        cpus: Number of threads for BWA.
        mincovpct: Coverage cutoff as a percent of the average N50-contig coverage.
        debug: Print per-contig details and keep the work directory.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.
    """
    workdir, custom_workdir = make_workdir(workdir, "covpurge")
    try:
        bamthreads = min(cpus, 4)
        forward_reads = str(Path(read1).resolve())
        reverse_reads = str(Path(read2).resolve()) if read2 else None

        assembly_working = "assembly.fasta"
        blob_bam = "remapped.bam"
        assembly_path = str(Path(workdir, assembly_working))
        shutil.copyfile(input, assembly_path)
        num_seqs, assembly_size = fasta_stats(assembly_path)
        logger.info(f"Assembly is {num_seqs:,} contigs and {assembly_size:,} bp")

        # reuse the BAM of a previous run in the same work directory
        if not Path(workdir, blob_bam).is_file():
            logger.info("Building BWA index")
            run_cmd(["bwa", "index", assembly_working], debug, cwd=workdir)
            bwa_cmd = ["bwa", "mem", "-t", str(cpus), assembly_working, forward_reads]
            if reverse_reads:
                bwa_cmd.append(reverse_reads)
            logger.info("Aligning reads to assembly with BWA")
            align_to_sorted_bam(bwa_cmd, str(Path(workdir, blob_bam)), bamthreads, cwd=workdir, debug=debug)

        logger.info("Calculating read coverage per contig")
        fasta_bed = str(Path(workdir, "assembly.bed"))
        lengths = []
        with open(fasta_bed, "w") as bedout, open(assembly_path) as seq_in:
            for record in SeqIO.parse(seq_in, "fasta"):
                bedout.write(f"{record.id}\t{0}\t{len(record.seq)}\n")
                lengths.append(len(record.seq))

        n50, _ = calc_nx(lengths)
        coverage = {}
        coverage_bed = str(Path(workdir, "coverage.bed"))
        cov_cmd = ["samtools", "bedcov", Path(fasta_bed).name, blob_bam]
        with open(coverage_bed, "w") as bed_out:
            for line in execute(cov_cmd, workdir, debug):
                bed_out.write(line)
                if not line or line.startswith("\n") or line.count("\t") < 3:
                    continue
                cols = line.strip().split("\t")
                coverage[cols[0]] = (int(cols[2]), int(cols[3]) / float(cols[2]))

        # average coverage of the N50 contigs is the "true genome coverage" baseline
        n50_cov = []
        for contig, (contig_len, cov) in coverage.items():
            if debug:
                logger.info(f"{contig}; Len: {contig_len}; Cov: {cov:.2f}")
            if contig_len >= n50:
                n50_cov.append(cov)
        contigs_to_drop = set()
        if not n50_cov:
            logger.warning(f"No contig of at least the N50 length ({n50:,} bp) has coverage data; skipping the low-coverage filter")
        else:
            n50_avg_cov = sum(n50_cov) / len(n50_cov)
            min_coverage = float(n50_avg_cov * mincovpct / 100)
            logger.info(f"Average coverage for N50 contigs is {int(n50_avg_cov)}X")
            contigs_to_drop = {contig for contig, (_, cov) in coverage.items() if cov <= min_coverage}
            logger.info(f"Found {len(contigs_to_drop):,} contigs with coverage less than {min_coverage:.2f}X ({mincovpct}%)")
        if debug:
            logger.info("Contigs dropped due to coverage: {:}".format(",".join(sorted(contigs_to_drop))))

        num_seqs, assembly_size = filter_fasta(assembly_path, outfile, lambda seq_id: seq_id not in contigs_to_drop)
        logger.info(f"Coverage-purged assembly is {num_seqs:,} contigs and {assembly_size:,} bp")
        cleanup_workdir(workdir, debug, custom_workdir)

        if not pipe:
            next_out = next_step_name(outfile, ".rmdup.fasta")
            logger.info(f"Your next command might be:\nAAFTF rmdup -i {outfile} -o {next_out} -c {cpus}")
    finally:
        close_step_log()
