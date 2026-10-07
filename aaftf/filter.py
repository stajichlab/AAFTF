"""Remove reads that match a contamination database before assembly.

The database combines PhiX and UniVec (from ``AAFTF database``), the Illumina sequencing-artifact
and lambda phage FASTAs bundled with BBTools, and any sequences the user adds. Reads are screened
with BBDuk k-mer matching, or aligned with bowtie2, bwa or minimap2.
"""

import logging
import shutil
from pathlib import Path
from typing import Any

from aaftf.resources import SEQ_DBS
from aaftf.utility import align_to_sorted_bam, basename_from_reads, cleanup_workdir, close_step_log, concat_files, count_fastq, db_file, download_file, make_workdir, require_databases, run_cmd

__all__ = ["run", "bbtools_screens"]


logger = logging.getLogger(__name__)

ALIGNERS = ("bbduk", "bowtie2", "bwa", "minimap2")

# Illumina sequencing artifacts and the lambda phage genome (a common spike-in), bundled with
# BBTools in its resources/ folder; added to every contamination database so that all
# --aligner choices screen the same sequences
BBTOOLS_SCREENS = ("sequencing_artifacts.fa.gz", "lambda.fa.gz")


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
    """Remove reads matching the contamination database.

    Builds ``contamdb.fa`` in the work directory from PhiX, UniVec, the BBTools artifact and lambda
    FASTAs and any extra accessions, URLs or local FASTA files. With ``aligner="bbduk"`` reads are
    k-mer filtered by BBDuk (paired reads are interleaved first and de-interleaved afterwards);
    otherwise reads are aligned with bowtie2, bwa or minimap2 and the unmapped reads are extracted
    with ``samtools fastq``. Cleaned reads are written to ``{basename}_filtered_1.fastq.gz`` (and
    ``_2`` for paired reads) in the current directory.

    In a reused work directory, ``contamdb.fa`` is rebuilt when its sources change, and the
    aligner index and BAM file are rebuilt when they are older than ``contamdb.fa``.

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
        RuntimeError: If the aligner produces no BAM file.
    """
    if aligner not in ALIGNERS:
        raise ValueError(f"Unknown aligner {aligner}; use {', '.join(ALIGNERS)}")
    if not read1:
        raise ValueError("Must provide -1/--read1 FASTQ reads")
    workdir, custom_workdir = make_workdir(workdir, "filter")
    try:
        forward_reads = str(Path(read1).resolve())
        reverse_reads = str(Path(read2).resolve()) if read2 else None
        basename = basename or basename_from_reads(forward_reads)
        clean_1 = f"{basename}_filtered_1.fastq.gz"
        clean_2 = f"{basename}_filtered_2.fastq.gz" if reverse_reads else None

        contamdb = _build_contamdb(workdir, screen_accessions, screen_urls, screen_local)

        pairs = 2 if reverse_reads else 1
        total = count_fastq(forward_reads) * pairs
        logger.info(f"Loading {total:,} total reads")

        if aligner == "bbduk":
            _filter_bbduk(contamdb, forward_reads, reverse_reads, clean_1, clean_2, workdir, cpus, memory, debug)
        else:
            align_bam = str(Path(workdir, basename + "_contam_db.bam"))
            _align(aligner, contamdb, forward_reads, reverse_reads, align_bam, workdir, cpus, debug)
            # keep only reads (pairs: both mates) that did not map to the contamination database
            if reverse_reads:
                samtools_cmd = ["samtools", "fastq", "-f", "12", "-1", clean_1, "-2", clean_2, align_bam]
            else:
                samtools_cmd = ["samtools", "fastq", "-f", "4", "-0", clean_1, align_bam]
            run_cmd(samtools_cmd, debug)
        cleanup_workdir(workdir, debug, custom_workdir)

        clean = count_fastq(clean_1) * pairs
        logger.info(f"{total - clean:,} reads matched the contamination database; {clean:,} reads kept")
        if reverse_reads:
            logger.info(f"Filtering complete:\nFor: {clean_1}\nRev: {clean_2}")
        else:
            logger.info(f"Filtering complete:\nSingle: {clean_1}")
        if not pipe:
            reads = f"-1 {clean_1} -2 {clean_2}" if clean_2 else f"-1 {clean_1}"
            logger.info(f"Your next command might be:\nAAFTF assemble --method spades {reads} -c {cpus} -o {basename}.spades.fasta")
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


def _build_contamdb(workdir: str, screen_accessions: list[str] | None, screen_urls: list[str] | None, screen_local: list[str] | None) -> str:
    """Write ``<workdir>/contamdb.fa`` from every screening source, unless it is already current.

    The sources are PhiX and UniVec (``AAFTF database``), the BBTools artifact and lambda FASTAs,
    then the extra accessions (downloaded once into the database folder), URLs and local files.
    The list of sources, with their sizes and modification times, is kept next to the database in
    ``contamdb.fa.sources``; the database is rebuilt whenever that list changes.

    Args:
        workdir: Work directory to write the database into.
        screen_accessions: GenBank nucleotide accessions to add.
        screen_urls: URLs of FASTA files to add.
        screen_local: Local FASTA files to add.

    Returns:
        Path of ``contamdb.fa``.
    """
    sources = require_databases(["phix", "univec"]) + bbtools_screens()
    sources += [download_file(SEQ_DBS["nucleotide"] % acc, db_file(acc + ".fna")) for acc in screen_accessions or []]
    sources += [download_file(url, str(Path(workdir, Path(url).name))) for url in screen_urls or []]
    sources += [str(Path(f).resolve()) for f in screen_local or []]

    contamdb = Path(workdir, "contamdb.fa")
    manifest = Path(workdir, "contamdb.fa.sources")
    listing = "".join(f"{src}\t{Path(src).stat().st_size}\t{Path(src).stat().st_mtime_ns}\n" for src in sources)
    if contamdb.is_file() and manifest.is_file() and manifest.read_text() == listing:
        logger.info(f"Reusing the contamination database {contamdb}")
    else:
        logger.info("Generating combined contamination database {} from:\n{}".format(contamdb, "\n".join(sources)))
        concat_files(sources, contamdb)
        manifest.write_text(listing)
    return str(contamdb)


def _filter_bbduk(contamdb: str, forward_reads: str, reverse_reads: str | None, clean_1: str, clean_2: str | None, workdir: str, cpus: int, memory: int, debug: bool) -> None:
    """K-mer filter the reads against ``contamdb`` with BBDuk into ``clean_1`` (and ``clean_2``).

    Paired reads are interleaved first and de-interleaved afterwards: BBDuk's paired mode
    (``in=``/``in2=``) silently truncates large paired FASTQ in this BBTools build, while its
    single-file reader handles an interleaved file correctly.
    """
    bbduk = ["bbduk.sh", f"-Xmx{memory}g", f"t={cpus}", "hdist=1", "k=27", "overwrite=true", f"ref={contamdb}"]
    logger.info("Kmer filtering reads using BBDuk")
    if reverse_reads:
        interleaved_in = str(Path(workdir, "reads_ivl.fq.gz"))
        interleaved_out = str(Path(workdir, "reads_ivl.clean.fq.gz"))
        run_cmd(["shuffle.sh", f"in1={forward_reads}", f"in2={reverse_reads}", f"out={interleaved_in}"], debug)
        run_cmd(bbduk + [f"in={interleaved_in}", "interleaved=true", f"out={interleaved_out}"], debug)
        run_cmd(["reformat.sh", f"in={interleaved_out}", f"out1={clean_1}", f"out2={clean_2}"], debug)
    else:
        run_cmd(bbduk + [f"in={forward_reads}", f"out={clean_1}"], debug)


def _align(aligner: str, contamdb: str, forward_reads: str, reverse_reads: str | None, align_bam: str, workdir: str, cpus: int, debug: bool) -> None:
    """Align the reads to ``contamdb`` into the sorted ``align_bam``, unless a current one exists.

    Raises:
        RuntimeError: If no BAM file was produced.
    """
    if _is_current(align_bam, contamdb):
        logger.info(f"Reusing the alignment {align_bam}")
        return
    contam = Path(contamdb).name  # the commands run in the work directory
    if aligner == "bowtie2":
        _build_index(contamdb + ".1.bt2", contamdb, ["bowtie2-build", contamdb, contamdb], debug)
        cmd = ["bowtie2", "-x", contam, "-p", str(cpus), "--very-sensitive"]
        cmd += ["-1", forward_reads, "-2", reverse_reads] if reverse_reads else ["-U", forward_reads]
    elif aligner == "bwa":
        _build_index(contamdb + ".amb", contamdb, ["bwa", "index", contamdb], debug)
        cmd = ["bwa", "mem", "-t", str(cpus), contam, forward_reads] + ([reverse_reads] if reverse_reads else [])
    else:
        cmd = ["minimap2", "-ax", "sr", "-t", str(cpus), contam, forward_reads] + ([reverse_reads] if reverse_reads else [])
    logger.info(f"Aligning reads to contamination database using {aligner}")
    align_to_sorted_bam(cmd, align_bam, min(cpus, 4), cwd=workdir, debug=debug)
    if not Path(align_bam).is_file():
        raise RuntimeError(f"{aligner} alignment produced no BAM file {align_bam} (rerun with -v to see the {aligner} output)")


def _build_index(marker_file: str, contamdb: str, build_cmd: list[str], debug: bool) -> None:
    """Run ``build_cmd`` unless ``marker_file`` (an index file it writes) is newer than ``contamdb``."""
    if not _is_current(marker_file, contamdb):
        run_cmd(build_cmd, debug, quiet_stdout=True)


def _is_current(target: str, source: str) -> bool:
    """Return whether ``target`` exists and was modified no earlier than ``source``."""
    return Path(target).is_file() and Path(target).stat().st_mtime >= Path(source).stat().st_mtime
