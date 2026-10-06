"""Run a genome assembly using defaults suited to fungi.

The assembler is chosen with --method (no default): SPAdes, MEGAHIT or Unicycler for
Illumina reads, or Flye for long reads. There is some access to updating parameters but
this entire package is intended to be a general solution for draft genome processing en masse.
"""

import logging
import os
import re
import shutil
import uuid
from pathlib import Path
from typing import Any

from aaftf.utility import fasta_stats, run_cmd, safe_remove

__all__ = ["run", "run_spades", "run_megahit", "run_unicycler", "run_flye"]

# assemble --longread_type -> Flye read-type flag
FLYE_READ_TYPES = {"nano-raw": "--nano-raw", "nano-hq": "--nano-hq", "pacbio-raw": "--pacbio-raw", "pacbio-hifi": "--pacbio-hifi"}


logger = logging.getLogger(__name__)


def run(
    out: str | None,
    method: str,
    read1: str | None = None,
    workdir: str | None = None,
    cpus: int = 1,
    memory: int = 32,
    isolate: bool = True,
    careful: bool = True,
    assembler_args: list[str] | None = None,
    tmpdir: str | None = None,
    read2: str | None = None,
    longreads: str | None = None,
    merged: str | None = None,
    longread_type: str = "nano-hq",
    genome_size: str | None = None,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the ``assemble`` subcommand by dispatching to the chosen assembler.

    Args:
        out: Output assembly FASTA; derived from ``read1`` (``longreads`` for Flye) if None.
        method: Assembler: ``"spades"``, ``"megahit"``, ``"unicycler"`` or ``"flye"``; there is no
            default (``"masurca"``/``"nextdenovo"`` only log that they are not implemented).
        read1: Read 1 (forward, or single-end) FASTQ; required by every method except ``flye``.
        workdir: Assembler output directory; a unique name is generated if None.
        cpus: Number of threads.
        memory: Max memory in GB (SPAdes ``--mem``; converted to bytes for MEGAHIT ``--memory``).
        isolate: Pass ``--isolate`` to SPAdes.
        careful: Pass ``--careful`` to SPAdes (only when ``isolate`` is False).
        assembler_args: Extra arguments appended to the assembler command.
        tmpdir: Assembler temporary directory.
        read2: Read 2 (reverse) FASTQ, or None for single-end.
        longreads: Long-read FASTQ (Unicycler, and required by Flye).
        merged: Merged-pair FASTQ, or None.
        longread_type: Flye read type (a ``FLYE_READ_TYPES`` key).
        genome_size: Estimated genome size for Flye (e.g. ``"40m"``), or None.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``); ignored.

    Raises:
        ValueError: If ``method`` is not a known assembler, or ``flye`` is given no ``longreads``.
    """
    if method == "flye":
        run_flye(workdir=workdir, cpus=cpus, longreads=longreads, longread_type=longread_type, genome_size=genome_size, assembler_args=assembler_args, out=out, debug=debug, pipe=pipe)
    elif method == "spades":
        run_spades(workdir=workdir, cpus=cpus, memory=memory, isolate=isolate, careful=careful, assembler_args=assembler_args, tmpdir=tmpdir, read1=read1, read2=read2, merged=merged, out=out, debug=debug, pipe=pipe)
    elif method == "megahit":
        run_megahit(workdir=workdir, cpus=cpus, memory=memory, assembler_args=assembler_args, tmpdir=tmpdir, read1=read1, read2=read2, merged=merged, out=out, debug=debug, pipe=pipe)
    elif method == "masurca":
        logger.info("Masurca assembly is not yet implemented in AAFTF")
    elif method == "nextdenovo":
        logger.info("NextDenovo assembly is not yet implemented in AAFTF")
    elif method == "unicycler":
        run_unicycler(workdir=workdir, cpus=cpus, read1=read1, read2=read2, longreads=longreads, merged=merged, out=out, debug=debug, pipe=pipe)
    else:
        raise ValueError(f"Unknown assembler method {method}")


def run_spades(
    workdir: str | None = None,
    cpus: int = 1,
    memory: int = 32,
    isolate: bool = True,
    careful: bool = True,
    assembler_args: list[str] | None = None,
    tmpdir: str | None = None,
    read1: str | None = None,
    read2: str | None = None,
    merged: str | None = None,
    out: str | None = None,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the SPAdes assembler and copy ``scaffolds.fasta`` to the output file.

    If ``workdir`` already exists, SPAdes is restarted from its last checkpoint instead.

    Args:
        workdir: SPAdes output directory; a unique ``spades_*`` name is generated if None (and removed after the run unless debug).
        cpus: Number of threads.
        memory: Memory limit in GB.
        isolate: Pass ``--isolate``.
        careful: Pass ``--careful`` (only when ``isolate`` is False).
        assembler_args: Extra SPAdes arguments.
        tmpdir: SPAdes temporary directory.
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ, or None.
        merged: Merged-pair FASTQ, or None.
        out: Output assembly FASTA; derived from ``read1`` if None.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Extra keyword arguments; ignored.
    """
    auto_workdir = not workdir
    if not workdir:
        workdir = "spades_" + str(uuid.uuid4())[:8]

    runcmd = ["spades.py", "--threads", str(cpus), "--mem", str(memory), "-o", workdir]

    if isolate:
        runcmd.extend(["--isolate"])
    elif careful:
        runcmd.extend(["--careful"])

    if assembler_args:
        runcmd.extend(assembler_args)

    if "--meta" not in runcmd:
        runcmd.extend(["--cov-cutoff", "auto"])

    if tmpdir:
        runcmd.extend(["--tmp-dir", tmpdir])

    forward_reads, reverse_reads = _resolve_reads(read1, read2)

    if not reverse_reads:
        runcmd.extend(["-s", forward_reads])
        if merged:
            # no paired library to attach merged reads to: add them as a second single-read library
            runcmd.extend(["--s", "2", merged])
    else:
        runcmd.extend(["-1", forward_reads, "-2", reverse_reads])
        if merged:
            # merged pairs belong to the paired-end library (same as --pe-m 1)
            runcmd.extend(["--merged", merged])

    # this basically overrides everything above and only runs --restart-from option
    if Path(workdir).is_dir():
        runcmd = ["spades.py", "-o", workdir, "--threads", str(cpus), "--mem", str(memory), "--restart-from", "last"]

    logger.info("Assembling FASTQ data using Spades")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, forward_reads, ".spades.fasta")
    _finish_assembly(Path(workdir, "scaffolds.fasta"), final_out, "Spades", cpus, pipe)
    if auto_workdir and not debug:
        safe_remove(workdir)


def run_megahit(
    workdir: str | None = None,
    cpus: int = 1,
    memory: int | None = None,
    assembler_args: list[str] | None = None,
    tmpdir: str | None = None,
    read1: str | None = None,
    read2: str | None = None,
    merged: str | None = None,
    out: str | None = None,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the MEGAHIT assembler, which is faster but may be less accurate than SPAdes.

    Args:
        workdir: MEGAHIT output directory; ``megahit_<pid>`` if None (and removed after the run unless debug).
        cpus: Number of threads.
        memory: Max memory in GB, or None to use the MEGAHIT default (90% of RAM).
        assembler_args: Extra MEGAHIT arguments.
        tmpdir: Temporary directory.
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ, or None.
        merged: Merged-pair FASTQ, or None. MEGAHIT has no merged-read option, so it is passed with
            ``-r`` as a single-end library.
        out: Output assembly FASTA; derived from ``read1`` if None.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Extra keyword arguments; ignored.
    """
    auto_workdir = not workdir
    if not workdir:
        workdir = "megahit_" + str(os.getpid())

    runcmd = ["megahit", "-t", str(cpus), "-o", workdir]

    if assembler_args:
        runcmd.extend(assembler_args)

    if memory:
        runcmd.extend(["--memory", str(memory * 10**9)])  # MEGAHIT takes bytes

    if tmpdir:
        runcmd.extend(["--tmp-dir", tmpdir])

    forward_reads, reverse_reads = _resolve_reads(read1, read2)

    single_reads = [forward_reads] if not reverse_reads else []
    if merged:
        single_reads.append(merged)
    if reverse_reads:
        runcmd.extend(["-1", forward_reads, "-2", reverse_reads])
    if single_reads:
        runcmd.extend(["-r", ",".join(single_reads)])

    if Path(workdir).is_dir():
        raise FileExistsError(f"MEGAHIT output folder {workdir} already exists (MEGAHIT cannot resume); remove it or pass another -w/--workdir")

    logger.info("Assembling FASTQ data using megahit")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, forward_reads, ".megahit.fasta")
    _finish_assembly(Path(workdir, "final.contigs.fa"), final_out, "Megahit", cpus, pipe)
    if auto_workdir and not debug:
        safe_remove(workdir)


def run_unicycler(
    workdir: str | None = None,
    cpus: int = 1,
    read1: str | None = None,
    read2: str | None = None,
    longreads: str | None = None,
    merged: str | None = None,
    out: str | None = None,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the Unicycler assembler.

    Args:
        workdir: Unicycler output directory; a unique ``unicycler_*`` name is generated if None (and removed after the run unless debug).
        cpus: Number of threads.
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ, or None.
        longreads: Long-read FASTQ passed as ``--long``, or None.
        merged: Merged-pair FASTQ, passed as ``--unpaired`` alongside paired reads (Unicycler takes a
            single ``--unpaired`` file, so it is ignored for single-end input).
        out: Output assembly FASTA; derived from ``read1`` if None.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Extra keyword arguments; ignored.
    """
    auto_workdir = not workdir
    if not workdir:
        workdir = "unicycler_" + str(uuid.uuid4())[:8]

    runcmd = ["unicycler", "--threads", str(cpus), "-o", workdir]

    # if memory:
    #    runcmd.extend(['--spades_options', f'-m {memory}'])

    forward_reads, reverse_reads = _resolve_reads(read1, read2)

    if longreads:
        runcmd.extend(["--long", longreads])

    if not reverse_reads:
        runcmd.extend(["--unpaired", forward_reads])
    else:
        runcmd.extend(["--short1", forward_reads, "--short2", reverse_reads])
        if merged:
            runcmd.extend(["--unpaired", merged])

    # not supporting restarting a run
    # this basically overrides everything above and only runs --restart-from option
    #    if Path(workdir).is_dir():
    #    runcmd = ['unicycler', '-o', workdir,
    #            '--threads', str(cpus),
    #            '--mem', memory,
    #            '--restart-from last']

    logger.info("Assembling FASTQ data using Unicycler")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, forward_reads, ".unicycler.fasta")
    _finish_assembly(Path(workdir, "assembly.fasta"), final_out, "Unicycler", cpus, pipe)
    if auto_workdir and not debug:
        safe_remove(workdir)


def run_flye(
    workdir: str | None = None,
    cpus: int = 1,
    longreads: str | None = None,
    longread_type: str = "nano-hq",
    genome_size: str | None = None,
    assembler_args: list[str] | None = None,
    out: str | None = None,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the Flye long-read assembler.

    If ``workdir`` holds a previous Flye run (its ``params.json``), Flye is restarted with ``--resume``.

    Args:
        workdir: Flye output directory; a unique ``flye_*`` name is generated if None (and removed after the run unless debug).
        cpus: Number of threads.
        longreads: Long-read FASTQ.
        longread_type: Read type, a ``FLYE_READ_TYPES`` key (e.g. ``"nano-hq"``).
        genome_size: Estimated genome size (e.g. ``"40m"``), or None to let Flye estimate it.
        assembler_args: Extra Flye arguments.
        out: Output assembly FASTA; derived from ``longreads`` if None.
        debug: Show external command output when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Extra keyword arguments; ignored.

    Raises:
        ValueError: If ``longreads`` is not given or ``longread_type`` is unknown.
    """
    if not longreads:
        raise ValueError("Flye needs long reads, provide -lr/--longreads")
    if longread_type not in FLYE_READ_TYPES:
        raise ValueError(f"Unknown long-read type {longread_type}; choose from {', '.join(FLYE_READ_TYPES)}")
    longreads = str(Path(longreads).resolve())
    auto_workdir = not workdir
    if not workdir:
        workdir = "flye_" + str(uuid.uuid4())[:8]

    runcmd = ["flye", FLYE_READ_TYPES[longread_type], longreads, "--out-dir", workdir, "--threads", str(cpus)]
    if genome_size:
        runcmd.extend(["--genome-size", genome_size])
    if assembler_args:
        runcmd.extend(assembler_args)
    if Path(workdir, "params.json").is_file():
        logger.info(f"Flye run found in {workdir}, resuming it")
        runcmd.append("--resume")

    logger.info("Assembling long reads using Flye")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, longreads, ".flye.fasta")
    _finish_assembly(Path(workdir, "assembly.fasta"), final_out, "Flye", cpus, pipe)
    if auto_workdir and not debug:
        safe_remove(workdir)


def _resolve_reads(read1: str | None, read2: str | None) -> tuple[str, str | None]:
    """Resolve absolute paths for the forward and reverse reads.

    Args:
        read1: Forward reads path.
        read2: Reverse reads path, or None.

    Returns:
        Tuple of (absolute forward path, absolute reverse path or None).

    Raises:
        ValueError: If ``read1`` is not given.
    """
    forward_reads = str(Path(read1).resolve()) if read1 else None
    reverse_reads = str(Path(read2).resolve()) if read2 else None
    if not forward_reads:
        raise ValueError("Unable to locate FASTQ raw reads, provide --read1")
    return forward_reads, reverse_reads


def _derive_final_out(out: str | None, forward_reads: str, suffix: str) -> str:
    """Derive the assembly output FASTA filename from ``out``, or from the input read filename.

    Args:
        out: User-supplied output path; returned unchanged if set.
        forward_reads: Forward reads path, whose name minus ``.fastq``/``.fq`` extension is the prefix.
        suffix: Suffix appended to the prefix (e.g. ``".spades.fasta"``).

    Returns:
        The output FASTA filename.
    """
    if out:
        return out
    prefix = Path(forward_reads).name
    m = re.search(r"(\S+)\.(fastq|fq)(\.\S+)?", prefix)
    if m:
        prefix = m.group(1)
    return prefix + suffix


def _finish_assembly(src: str | Path, final_out: str, tool_name: str, cpus: int, pipe: bool) -> None:
    """Copy the assembler's raw output to ``final_out``, report stats, and log the next-step hint.

    Args:
        src: Assembler output FASTA.
        final_out: Destination FASTA path.
        tool_name: Assembler name used in log messages.
        cpus: Thread count shown in the suggested command.
        pipe: Suppress the "next command" hint.

    Raises:
        RuntimeError: If the assembler did not produce ``src``.
    """
    if not Path(src).is_file():
        raise RuntimeError(f"{tool_name} assembly output {src} is missing -- check the {tool_name} log in {Path(src).parent}")
    shutil.copyfile(str(src), final_out)
    logger.info(f"{tool_name} assembly finished: {final_out}")
    num_seqs, assembly_size = fasta_stats(final_out)
    logger.info(f"Assembly is {num_seqs:,} scaffolds and {assembly_size:,} bp")

    if not pipe:
        logger.info(f"Your next command might be:\nAAFTF vecscreen -i {final_out} -c {cpus}")
