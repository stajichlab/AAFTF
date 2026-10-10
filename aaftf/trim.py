"""Adapter- and quality-trim FASTQ reads.

This is usually for Illumina reads. Trimming uses BBDuk (default), fastp
(optionally merging paired reads) or Trimmomatic; Trimmomatic expects the
adaptor sequence files to be in its installed folder.
"""

import logging
import shutil
from pathlib import Path
from typing import Any

from aaftf.utility import basename_from_reads, count_fastq, run_cmd, safe_remove

__all__ = ["TRIMMOMATIC_TRUSEQSE", "TRIMMOMATIC_TRUSEQPE", "run", "run_bbduk", "run_trimmomatic", "run_fastp"]


TRIMMOMATIC_TRUSEQSE = "adapters/TruSeq3-SE.fa"

TRIMMOMATIC_TRUSEQPE = "adapters/TruSeq3-PE.fa"

# process trimming reads with trimmomatic
# Homebrew install of trimmomatic uses a shell script
"""
#!/bin/bash -l
TRIMJAR=/usr/local/Cellar/trimmomatic/0.36/libexec/trimmomatic-0.36.jar
exec java -jar $TRIMJAR "$@"
"""

# while bioconda install uses a python script that launches java apps

logger = logging.getLogger(__name__)


# flake8: noqa: C901
def run(
    read1: str,
    read2: str | None = None,
    basename: str | None = None,
    method: str = "bbduk",
    cpus: int = 1,
    memory: int = 8,
    minlen: int = 75,
    avgqual: int = 10,
    trimmomatic_adaptors: str = "TruSeq3-PE.fa",
    trimmomatic_clip: str = "2:30:10",
    trimmomatic_leadingwindow: int = 3,
    trimmomatic_trailingwindow: int = 3,
    trimmomatic_slidingwindow: str = "4:15",
    trimmomatic_quality: str = "phred33",
    merge: bool = False,
    dedup: bool = False,
    cutfront: bool = False,
    cuttail: bool = False,
    cutright: bool = False,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the ``trim`` subcommand: count input reads and dispatch to the chosen trimmer.

    Args:
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ for paired-end data, or None.
        basename: Output file prefix; derived from ``read1`` if not given.
        method: Trimmer to use: ``"bbduk"``, ``"trimmomatic"`` or ``"fastp"``.
        cpus: Number of threads.
        memory: Java heap size in GB (BBDuk only).
        minlen: Minimum read length after trimming.
        avgqual: Minimum average read quality (BBDuk and fastp).
        trimmomatic_adaptors: Trimmomatic adaptor FASTA.
        trimmomatic_clip: Trimmomatic ILLUMINACLIP settings.
        trimmomatic_leadingwindow: Trimmomatic LEADING quality.
        trimmomatic_trailingwindow: Trimmomatic TRAILING quality.
        trimmomatic_slidingwindow: Trimmomatic SLIDINGWINDOW settings.
        trimmomatic_quality: Quality encoding, ``"phred33"`` or ``"phred64"``.
        merge: Merge overlapping pairs (fastp only).
        dedup: Remove duplicate reads (fastp only).
        cutfront: Enable fastp ``--cut_front``.
        cuttail: Enable fastp ``--cut_tail``.
        cutright: Enable fastp ``--cut_right``.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``); ignored.
    """
    if not basename:
        basename = basename_from_reads(read1)

    total = count_fastq(read1)
    if read2:
        total = total * 2
    logger.info(f"Loading {total:,} total reads")

    if method == "bbduk":
        run_bbduk(read1, read2, basename, cpus, memory, minlen, avgqual, debug, pipe)
    elif method == "trimmomatic":
        run_trimmomatic(
            read1,
            read2,
            basename,
            cpus,
            minlen,
            trimmomatic_adaptors,
            trimmomatic_clip,
            trimmomatic_leadingwindow,
            trimmomatic_trailingwindow,
            trimmomatic_slidingwindow,
            trimmomatic_quality,
            debug,
            pipe,
        )
    elif method == "fastp":
        run_fastp(read1, read2, basename, cpus, minlen, avgqual, merge, dedup, cutfront, cuttail, cutright, debug, pipe)
    else:
        logger.info(f"Unknown trimming method: {method}")


def run_bbduk(read1: str, read2: str | None, basename: str, cpus: int, memory: int, minlen: int, avgqual: int, debug: bool, pipe: bool) -> None:
    """Trim reads with BBDuk.

    Paired reads are interleaved with ``shuffle.sh`` before trimming and split back into
    ``<basename>_1P``/``_2P`` files with ``reformat.sh``; single-end output is ``<basename>_1U``.

    Args:
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ, or None for single-end.
        basename: Output file prefix.
        cpus: Number of threads.
        memory: Java heap size in GB.
        minlen: Minimum read length after trimming.
        avgqual: Minimum average read quality (``maq``).
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint.
    """
    java_mem = f"-Xmx{memory}g"

    logger.info("Adapter trimming using BBDuk")
    bbduk_base = [
        "bbduk.sh",
        java_mem,
        "ref=adapters",
        f"t={cpus}",
        "ktrim=r",
        "k=23",
        "mink=11",
        f"minlen={minlen}",
        "hdist=1",
        f"maq={avgqual}",
        "ftm=5",
        "tpe",
        "tbo",
        "overwrite=true",
    ]
    if read1 and read2:
        # Paired mode (in1=/in2=) hits a bug in this BBDuk build's
        # PairStreamer on large/variable-length paired FASTQ: it silently
        # truncates the stream after a few hundred reads instead of
        # erroring loudly. Feed it an INTERLEAVED single file instead
        # (bbduk's single-end reader handles the full file correctly),
        # then de-interleave the trimmed output.
        interleaved_in = f"{basename}_ivl.fq.gz"
        interleaved_out = f"{basename}_ivl.trimmed.fq.gz"
        shuffle_cmd = ["shuffle.sh", f"in1={read1}", f"in2={read2}", f"out={interleaved_in}"]
        run_cmd(shuffle_cmd, debug)

        cmd = bbduk_base + [f"in={interleaved_in}", "interleaved=true", f"out={interleaved_out}"]
        run_cmd(cmd, debug)

        reformat_cmd = [
            "reformat.sh",
            f"in={interleaved_out}",
            f"out1={basename}_1P.fastq.gz",
            f"out2={basename}_2P.fastq.gz",
        ]
        run_cmd(reformat_cmd, debug)
        safe_remove(interleaved_in)
        safe_remove(interleaved_out)
    elif read1:
        cmd = bbduk_base + [f"in={read1}", f"out={basename}_1U.fastq.gz"]
        run_cmd(cmd, debug)

    _report_trimmed(basename, read2, cpus, pipe)


def run_trimmomatic(
    read1: str,
    read2: str | None,
    basename: str,
    cpus: int,
    minlen: int,
    trimmomatic_adaptors: str,
    trimmomatic_clip: str,
    trimmomatic_leadingwindow: int,
    trimmomatic_trailingwindow: int,
    trimmomatic_slidingwindow: str,
    trimmomatic_quality: str,
    debug: bool,
    pipe: bool,
) -> None:
    """Trim reads with Trimmomatic.

    If ``trimmomatic_adaptors`` does not exist, the TruSeq3 adaptor file is searched for
    next to the jar and under parent ``share/trimmomatic`` directories; the function logs
    and returns without trimming if none is found.

    Args:
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ, or None for single-end.
        basename: Output file prefix.
        cpus: Number of threads.
        minlen: Minimum read length after trimming.
        trimmomatic_adaptors: Adaptor FASTA for ILLUMINACLIP.
        trimmomatic_clip: ILLUMINACLIP settings.
        trimmomatic_leadingwindow: LEADING quality.
        trimmomatic_trailingwindow: TRAILING quality.
        trimmomatic_slidingwindow: SLIDINGWINDOW settings.
        trimmomatic_quality: Quality encoding, ``"phred33"`` or ``"phred64"``.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint.

    Raises:
        FileNotFoundError: If the Trimmomatic jar or the adaptors file cannot be located.
        ValueError: If ``read1`` is not given.
    """
    if not read1:
        raise ValueError("Must provide read1 and read2 pairs or a single read set")
    trimmomatic_path = _find_trimmomatic()
    if trimmomatic_path:
        jarfile = trimmomatic_path
    else:
        raise FileNotFoundError("Trimmomatic cannot be found - please provide location of trimmomatic.jar file.")

    path_to_adaptors = trimmomatic_adaptors
    leadingwindow = f"LEADING:{trimmomatic_leadingwindow}"
    trailingwindow = f"TRAILING:{trimmomatic_trailingwindow}"
    slidingwindow = f"SLIDINGWINDOW:{trimmomatic_slidingwindow}"

    quality = trimmomatic_quality
    quality = f"-{quality}"  # add leading dash

    if not Path(path_to_adaptors).exists():
        adaptor_name = TRIMMOMATIC_TRUSEQPE if read2 else TRIMMOMATIC_TRUSEQSE
        path_to_adaptors = str(Path(jarfile).parent / adaptor_name)

        # otherwise look for <prefix>/share/trimmomatic/<adaptors> in each folder above the jar
        findpath = Path(jarfile).parent
        while not Path(path_to_adaptors).exists():
            path_to_adaptors = str(findpath / "share" / "trimmomatic" / adaptor_name)
            if findpath.parent == findpath:  # reached filesystem root
                break
            findpath = findpath.parent

        if not Path(path_to_adaptors).exists():
            raise FileNotFoundError(f"Cannot find the Trimmomatic adaptors file {trimmomatic_adaptors}; pass its path with --trimmomatic_adaptors")
    clipstr = f"ILLUMINACLIP:{path_to_adaptors}:{trimmomatic_clip}"

    if read1 and read2:
        cmd = [
            "java",
            "-jar",
            jarfile,
            "PE",
            "-threads",
            str(cpus),
            quality,
            read1,
            read2,
            basename + "_1P.fastq.gz",
            basename + "_1U.fastq.gz",
            basename + "_2P.fastq.gz",
            basename + "_2U.fastq.gz",
            clipstr,
            leadingwindow,
            trailingwindow,
            slidingwindow,
            f"MINLEN:{minlen}",
        ]
    else:
        cmd = [
            "java",
            "-jar",
            jarfile,
            "SE",
            "-threads",
            str(cpus),
            quality,
            read1,
            basename + "_1U.fastq.gz",
            clipstr,
            leadingwindow,
            trailingwindow,
            slidingwindow,
            f"MINLEN:{minlen}",
        ]

    logger.info("Running trimmomatic adapter and quality trimming")
    run_cmd(cmd, debug)
    if read2:
        safe_remove(basename + "_1U.fastq.gz")
        safe_remove(basename + "_2U.fastq.gz")
        logger.info("Trimming finished:\nFor: {:}\nRev {:}".format(basename + "_1P.fastq.gz", basename + "_2P.fastq.gz"))
        if not pipe:
            logger.info("Your next command might be:\n" + "AAFTF filter -1 {:} -2 {:} -o {:} -c {:}".format(basename + "_1P.fastq.gz", basename + "_2P.fastq.gz", basename, cpus))
    else:
        logger.info("Trimming finished:\nSingle: {:}".format(basename + "_1U.fastq.gz"))
        if not pipe:
            logger.info("Your next command might be:\n" + "AAFTF filter -1 {:} -o {:} -c {:}".format(basename + "_1U.fastq.gz", basename, cpus))


def run_fastp(read1: str, read2: str | None, basename: str, cpus: int, minlen: int, avgqual: int, merge: bool, dedup: bool, cutfront: bool, cuttail: bool, cutright: bool, debug: bool, pipe: bool) -> None:
    """Trim reads with fastp, writing HTML and JSON reports alongside the trimmed FASTQ.

    Args:
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ, or None for single-end.
        basename: Output file prefix.
        cpus: Number of threads.
        minlen: Minimum read length after trimming.
        avgqual: Minimum average read quality.
        merge: Merge overlapping pairs into ``<basename>_MG.fastq.gz`` (paired only).
        dedup: Remove duplicate reads.
        cutfront: Enable ``--cut_front``.
        cuttail: Enable ``--cut_tail``.
        cutright: Enable ``--cut_right``.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint.
    """
    logger.info("Adapter trimming using fastp")
    cmd = [
        "fastp",
        "--low_complexity_filter",
        "-l",
        f"{minlen}",
        "--average_qual",
        f"{avgqual}",
        "-w",
        f"{cpus}",
    ]

    if read1 and read2:
        # could add merging ...
        cmd += [
            f"--in1={read1}",
            f"--in2={read2}",
            f"--out1={basename}_1P.fastq.gz",
            f"--out2={basename}_2P.fastq.gz",
        ]
        if merge:
            cmd += ["--merge", f"--merged_out={basename}_MG.fastq.gz"]

    elif read1:
        cmd += [f"--in={read1}", f"--out={basename}_1U.fastq.gz"]
    if dedup:
        cmd += ["--dedup"]
    if cutfront:
        cmd += ["--cut_front"]
    if cuttail:
        cmd += ["--cut_tail"]
    if cutright:
        cmd += ["--cut_right"]

    cmd += [f"--html={basename}.fastp.html", f"--json={basename}.fastp.json"]
    run_cmd(cmd, debug)

    _report_trimmed(basename, read2, cpus, pipe)


def _report_trimmed(basename: str, read2: str | None, cpus: int, pipe: bool) -> None:
    """Log the number of reads left after trimming and the suggested next command.

    Args:
        basename: Output file prefix of the trimmed reads.
        read2: Read 2 (reverse) FASTQ; truthy means paired-end output is counted.
        cpus: Thread count shown in the suggested command.
        pipe: Suppress the "next command" hint.
    """
    if read2:
        clean = count_fastq(f"{basename}_1P.fastq.gz")
        clean = clean * 2
        logger.info(f"{clean:,} reads remaining and writing to file")
        logger.info("Trimming finished:\nFor: {:}\nRev {:}".format(basename + "_1P.fastq.gz", basename + "_2P.fastq.gz"))
        if not pipe:
            logger.info("Your next command might be:\n" + "AAFTF filter -1 {:} -2 {:} -o {:} -c {:}".format(basename + "_1P.fastq.gz", basename + "_2P.fastq.gz", basename, cpus))
    else:
        clean = count_fastq(f"{basename}_1U.fastq.gz")
        logger.info(f"{clean:,} reads remaining and writing to file")
        logger.info("Trimming finished:\nSingle: {:}".format(basename + "_1U.fastq.gz"))
        if not pipe:
            logger.info("Your next command might be:\n" + "AAFTF filter -1 {:} -o {:} -c {:}".format(basename + "_1U.fastq.gz", basename, cpus))


def _find_trimmomatic() -> str | None:
    """Locate the Trimmomatic jar file from the ``trimmomatic`` launcher on PATH.

    Handles the Homebrew bash wrapper (reads the jar from its ``exec java`` line) and the
    bioconda Python wrapper (``trimmomatic.jar`` next to the script).

    Returns:
        The jar path, or None if the launcher is missing, of an unknown type, or a bash wrapper
        without a ``.jar`` on an ``exec java`` line.
    """
    trim_path = shutil.which("trimmomatic")
    if trim_path:
        with open(str(Path(trim_path).resolve())) as trim_shell:
            first_line = trim_shell.readline()
            if "#!/bin/bash" in first_line:  # homebrew get jar location
                for line in trim_shell:
                    if line.startswith("exec java"):
                        items = line.split(" ")
                        for x in items:
                            if x.endswith(".jar"):
                                return x
            elif "#!/usr/bin/env python" in first_line:
                trimjardir = Path(trim_path).resolve().parent
                return str(trimjardir / "trimmomatic.jar")
    return None
