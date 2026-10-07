"""Filter runs routines to remove sequence reads matching contaminant DB.

This will match contaminant database and PhiX using read mapping kmer
tools. See resources.py for these defaults.
"""

import logging
import shutil
from pathlib import Path
from typing import Any

from aaftf.resources import SEQ_DBS
from aaftf.utility import align_to_sorted_bam, bam_read_count, basename_from_reads, cleanup_workdir, close_step_log, concat_files, count_fastq, db_file, download_file, make_workdir, require_databases, run_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)

# Illumina sequencing artifacts and the lambda phage genome (a common spike-in), bundled with
# BBTools in its resources/ folder; added to every contamination database so that all
# --aligner choices screen the same sequences
BBTOOLS_SCREENS = ("sequencing_artifacts.fa.gz", "lambda.fa.gz")


# flake8: noqa: C901
def run(
    read1: str,
    read2: str | None = None,
    workdir: str | None = None,
    cpus: int = 1,
    screen_accessions: list[str] | None = None,
    screen_urls: list[str] | None = None,
    screen_local: list[str] | None = None,
    basename: str | None = None,
    aligner: str = "bbduk",
    memory: int = 8,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Remove reads matching the contamination database (PhiX, UniVec and any extra sequences).

    Builds ``contamdb.fa`` in the work directory from the PhiX/UniVec databases plus any extra
    accessions, URLs or local FASTA files. With ``aligner="bbduk"`` reads are k-mer filtered by
    BBDuk (paired reads are interleaved first and de-interleaved afterwards); otherwise reads are
    aligned with bowtie2, bwa or minimap2 and the unmapped reads are extracted with
    ``samtools fastq``. Cleaned reads are written to ``{basename}_filtered_*.fastq.gz`` in the
    current directory.

    Args:
        read1: Forward (or single-end) FASTQ file.
        read2: Reverse FASTQ file for paired-end data.
        workdir: Working directory; a temporary one is created when None.
        cpus: Number of threads.
        screen_accessions: GenBank nucleotide accessions to download and add to the screen.
        screen_urls: URLs of FASTA files to download and add to the screen.
        screen_local: Local FASTA files to add to the screen.
        basename: Output file prefix; derived from ``read1`` when None.
        aligner: One of ``bbduk``, ``bowtie2``, ``bwa`` or ``minimap2``.
        memory: Java heap size in GB for BBDuk.
        debug: Show external tool output and keep the work directory.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.

    Raises:
        ValueError: If ``read1`` is not given or ``aligner`` is unknown.
    """
    if aligner not in ("bbduk", "bowtie2", "bwa", "minimap2"):
        raise ValueError(f"Unknown aligner {aligner}; use bbduk, bowtie2, bwa or minimap2")
    workdir, custom_workdir = make_workdir(workdir, "filter")
    try:
        bamthreads = min(cpus, 4)

        # PhiX and UniVec come from `AAFTF database`, artifacts and lambda from BBTools;
        # extra accessions/URLs are fetched here
        contam_filenames = require_databases(["phix", "univec"]) + bbtools_screens()

        for acc in screen_accessions or []:
            contam_filenames.append(download_file(SEQ_DBS["nucleotide"] % acc, db_file(acc + ".fna")))

        for url in screen_urls or []:
            contam_filenames.append(download_file(url, str(Path(workdir, Path(url).name))))

        contam_filenames.extend(str(Path(f).resolve()) for f in screen_local or [])

        contamdb = str(Path(workdir, "contamdb.fa"))
        filelist = "\n".join(contam_filenames)
        logger.info(f"Generating combined contamination database {contamdb} from:\n{filelist}")
        newest_source = max(Path(f).stat().st_ctime for f in contam_filenames)
        if not Path(contamdb).exists() or Path(contamdb).stat().st_ctime < newest_source:
            concat_files(contam_filenames, contamdb)

        # find reads
        forward_reads, reverse_reads = (None,) * 2
        if read1:
            forward_reads = str(Path(read1).resolve())
        if read2:
            reverse_reads = str(Path(read2).resolve())
        if not forward_reads:
            raise ValueError("Must provide --read1, unable to locate FASTQ reads")
        total = count_fastq(forward_reads)
        if reverse_reads:
            total = total * 2
        logger.info(f"Loading {total:,} total reads")

        # seems like this needs to be stripping trailing extension?
        if not basename:
            basename = basename_from_reads(forward_reads)

        # logger.info('Loading {:,} FASTQ reads'.format(count_fastq(forReads)))

        align_bam = str(Path(workdir, basename + "_contam_db.bam"))
        clean_reads = basename + "_filtered"
        if aligner == "bbduk":
            logger.info("Kmer filtering reads using BBDuk")
            java_mem = f"-Xmx{memory}g"
            leftcleanfname = f"{clean_reads}_1.fastq.gz"
            if reverse_reads:
                # Paired mode (in=/in2=) hits a bug in this BBDuk build's
                # PairStreamer on large paired FASTQ: it silently truncates the
                # stream after a few hundred reads (with or without threading),
                # instead of erroring loudly. Feed it an INTERLEAVED single file
                # instead (bbduk's single-end reader handles the full file
                # correctly), then de-interleave the cleaned output.
                interleaved_in = str(Path(workdir, f"{basename}_ivl.fq.gz"))
                interleaved_out = str(Path(workdir, f"{basename}_ivl.clean.fq.gz"))
                shuffle_cmd = ["shuffle.sh", f"in1={forward_reads}", f"in2={reverse_reads}", f"out={interleaved_in}"]
                run_cmd(shuffle_cmd, debug)

                cmd = ["bbduk.sh", java_mem, f"t={cpus}", "hdist=1", "k=27", "overwrite=true", f"in={interleaved_in}", "interleaved=true", f"out={interleaved_out}"]
                cmd.append(f"ref={contamdb}")
                run_cmd(cmd, debug)

                reformat_cmd = ["reformat.sh", f"in={interleaved_out}", f"out1={clean_reads}_1.fastq.gz", f"out2={clean_reads}_2.fastq.gz"]
                run_cmd(reformat_cmd, debug)
            else:
                cmd = ["bbduk.sh", java_mem, f"t={cpus}", "hdist=1", "k=27", "overwrite=true"]
                cmd.extend([f"in={forward_reads}", f"out={clean_reads}_U.fastq.gz"])
                leftcleanfname = f"{clean_reads}_U.fastq.gz"
                cmd.append(f"ref={contamdb}")
                # cmd.extend(['prealloc','qhdist=1'])
                run_cmd(cmd, debug)

            cleanup_workdir(workdir, debug, custom_workdir)

            clean = count_fastq(leftcleanfname)
            if reverse_reads:
                clean = clean * 2  # might want to actually count - but should be always 2x
            logger.info(f"{(total - clean):,} reads mapped to contamination database")
            logger.info(f"{clean:,} reads unmapped and writing to file")

            if reverse_reads:
                logger.info(f"Filtering complete:\nFor: {clean_reads}_1.fastq.gz\nRev: {clean_reads}_2.fastq.gz")
                if not pipe:
                    logger.info(f"Your next command might be:\nAAFTF assemble --method spades -1 {clean_reads}_1.fastq.gz -2 {clean_reads}_2.fastq.gz -c {cpus} -o {basename}.spades.fasta")

            else:
                logger.info(f"Filtering complete:\nSingle: {clean_reads}_U.fastq.gz")
                if not pipe:
                    logger.info(f"Your next command might be:\nAAFTF assemble --method spades -1 {clean_reads}_U.fastq.gz -c {cpus} -o {basename}.spades.fasta")

            return

        elif aligner == "bowtie2":
            # likely not used and less accurate than bbmap?
            if not Path(align_bam).is_file():
                logger.info("Aligning reads to contamination database using bowtie2")
                _rebuild_index_if_stale(contamdb + ".1.bt2", contamdb, ["bowtie2-build", contamdb, contamdb], debug)

                bowtie_cmd = ["bowtie2", "-x", Path(contamdb).name, "-p", str(cpus), "--very-sensitive"]
                if forward_reads and reverse_reads:
                    bowtie_cmd = bowtie_cmd + ["-1", forward_reads, "-2", reverse_reads]
                elif forward_reads:
                    bowtie_cmd = bowtie_cmd + ["-U", forward_reads]

                align_to_sorted_bam(bowtie_cmd, align_bam, bamthreads, cwd=workdir, debug=debug)

        elif aligner == "bwa":
            # likely less accurate than bbduk so may not be used
            if not Path(align_bam).is_file():
                logger.info("Aligning reads to contamination database using BWA")
                _rebuild_index_if_stale(contamdb + ".amb", contamdb, ["bwa", "index", contamdb], debug)

                bwa_cmd = ["bwa", "mem", "-t", str(cpus), Path(contamdb).name, forward_reads]
                if reverse_reads:
                    bwa_cmd.append(reverse_reads)

                align_to_sorted_bam(bwa_cmd, align_bam, bamthreads, cwd=workdir, debug=debug)

        elif aligner == "minimap2":
            # likely not used but may be useful for pacbio/nanopore?
            if not Path(align_bam).is_file():
                logger.info("Aligning reads to contamination database using minimap2")

                minimap2_cmd = ["minimap2", "-ax", "sr", "-t", str(cpus), Path(contamdb).name, forward_reads]
                if reverse_reads:
                    minimap2_cmd.append(reverse_reads)

                align_to_sorted_bam(minimap2_cmd, align_bam, bamthreads, cwd=workdir, debug=debug)

        if Path(align_bam).is_file():
            # display mapping stats in terminal
            mapped, unmapped = bam_read_count(align_bam)
            logger.info(f"{mapped:,} reads mapped to contamination database")
            logger.info(f"{unmapped:,} reads unmapped and writing to file")
            # now output unmapped reads from bamfile
            # this needs to be -f 5 so unmapped-pairs
            if forward_reads and reverse_reads:
                samtools_cmd = ["samtools", "fastq", "-f", "12", "-1", clean_reads + "_1.fastq.gz", "-2", clean_reads + "_2.fastq.gz", align_bam]
            elif forward_reads:
                samtools_cmd = ["samtools", "fastq", "-f", "4", "-0", clean_reads + ".fastq.gz", align_bam]
            run_cmd(samtools_cmd, debug)
            cleanup_workdir(workdir, debug, custom_workdir)

            if reverse_reads:
                logger.info(f"Filtering complete:\nFor: {clean_reads}_1.fastq.gz\nRev: {clean_reads}_2.fastq.gz")
                if not pipe:
                    logger.info(f"Your next command might be:\nAAFTF assemble --method spades -1 {clean_reads}_1.fastq.gz -2 {clean_reads}_2.fastq.gz -c {cpus} -o {basename}.spades.fasta")
            else:
                logger.info(f"Filtering complete:\nSingle: {clean_reads}.fastq.gz")
                if not pipe:
                    logger.info(f"Your next command might be:\nAAFTF assemble --method spades -1 {clean_reads}.fastq.gz -c {cpus} -o {basename}.spades.fasta")
        else:
            raise RuntimeError(f"{aligner} alignment produced no BAM file {align_bam} (rerun with -v to see the {aligner} output)")
    finally:
        close_step_log()


def bbtools_screens() -> list[str]:
    """Return the BBTools-bundled screening FASTAs (``BBTOOLS_SCREENS``) from its ``resources/`` folder.

    The folder sits next to the real ``bbduk.sh`` script (after resolving the symlink conda/pixi
    put on ``PATH``).

    Returns:
        Paths of the sequencing-artifact and lambda phage FASTAs.

    Raises:
        FileNotFoundError: If ``bbduk.sh`` is not on ``PATH`` or the files are not in its resources folder.
    """
    bbduk = shutil.which("bbduk.sh")
    if not bbduk:
        raise FileNotFoundError("bbduk.sh not found on PATH; filter needs BBTools (bbmap) for its bundled artifact and lambda phage sequences")
    resources = Path(bbduk).resolve().parent / "resources"
    paths = [resources / name for name in BBTOOLS_SCREENS]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"BBTools screening sequence(s) not found: {', '.join(missing)}")
    return [str(path) for path in paths]


def _rebuild_index_if_stale(marker_file: str, contamdb: str, build_cmd: list[str], debug: bool) -> None:
    """Run ``build_cmd`` if the index marker file is missing or older than ``contamdb``.

    Args:
        marker_file: An index file whose presence/ctime marks a built index.
        contamdb: The contamination FASTA the index is built from.
        build_cmd: Index-building command to run.
        debug: Passed to ``run_cmd`` to show tool output.
    """
    marker = Path(marker_file)
    if not marker.exists() or marker.stat().st_ctime < Path(contamdb).stat().st_ctime:
        run_cmd(build_cmd, debug, quiet_stdout=True)
