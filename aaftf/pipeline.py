"""Support pipelining of AAFTF to simplify all-in-one runs.

Three pipelines, one per kind of sequencing data:

* ``pipeline_short`` (Illumina): trim → filter → assemble → vecscreen → [taxpurge] → covpurge → rmdup →
  sort → assess.
* ``pipeline_long`` (ONT/PacBio): assemble (Flye) → vecscreen → [taxpurge] → rmdup → polish (Racon) →
  sort → assess (no covpurge: it needs Illumina reads).
* ``pipeline_hybrid`` (Illumina + long reads): trim → filter → assemble (Flye on the long reads) →
  vecscreen → [taxpurge] → covpurge → rmdup → polish (Racon, then Polypolish or pypolca with the Illumina reads)
  → sort → assess. With ``--method unicycler`` both read types go to Unicycler, which polishes its own
  assembly, so the polish steps are skipped.

Optionally, NCBI FCS-adaptor (``fcs_screen``, ``--fcs_screen``) also screens the contigs for
adaptors right after vecscreen. NCBI FCS-GX (``fcs_gx_purge``) is not a pipeline step: its
~470 GB database must be loaded into RAM, so it is run on its own on the finished assemblies.

All three share the same order: the contigs are cleaned (vecscreen, taxpurge, covpurge, rmdup) before
any polishing, so no time is spent polishing contigs that are then discarded. ``taxpurge`` runs only
when ``--phylum`` is given, so the pipelines also work when the organism's taxonomy is unknown.

Every step runs with its own ``AAFTF <step>`` defaults; only the options given to the pipeline (and
the file names that chain the steps together) override them.
"""

import argparse as ap
import functools
import logging
import shutil
from pathlib import Path
from types import ModuleType
from typing import Any

import aaftf.assemble as assemble
import aaftf.assess as assess
import aaftf.covpurge as covpurge
import aaftf.depth as depth
import aaftf.fcs_screen as fcs_screen
import aaftf.filter as aaftf_filter
import aaftf.polish as polish
import aaftf.rmdup as rmdup
import aaftf.sort as aaftf_sort
import aaftf.taxpurge as taxpurge
import aaftf.trim as trim
import aaftf.vecscreen as vecscreen
from aaftf.utility import check_file

# minimap2 preset (depth --longread_preset) for each Flye read type (--longread_type)
LONGREAD_PRESETS = {"nano-raw": "map-ont", "nano-hq": "map-ont", "pacbio-raw": "map-pb", "pacbio-hifi": "map-hifi"}

__all__ = ["run_short", "run_long", "run_hybrid"]


logger = logging.getLogger(__name__)


def run_short(
    read1: str,
    basename: str,
    phylum: list[str] | None = None,
    read2: str | None = None,
    cpus: int = 1,
    tmpdir: str | None = None,
    assembler_args: list[str] | None = None,
    method: str = "spades",
    memory: int | None = None,
    minlen: int = 75,
    screen_accessions: list[str] | None = None,
    screen_urls: list[str] | None = None,
    mincontiglen: int = 500,
    workdir: str | None = None,
    sourdb: list[str] | None = None,
    fcs_screen_run: bool = False,
    container_engine: str | None = None,
    mincovpct: int = 5,
    debug: bool = False,
    quiet: bool = False,
    **kwargs: Any,
) -> None:
    """Run the short-read (Illumina) pipeline, skipping steps whose output file already exists.

    Args:
        read1: Read 1 (forward, or single-end) raw FASTQ reads.
        basename: Prefix for every output file (``{basename}_1P.fastq.gz``, ``{basename}.final.fasta``, ...).
        phylum: Phyla whose sourmash matches are kept by ``taxpurge``; None skips ``taxpurge``.
        read2: Read 2 (reverse) raw FASTQ reads, or None for single-end data.
        cpus: Threads for every step that takes ``-c/--cpus``.
        tmpdir: Assembler temporary directory.
        assembler_args: Extra arguments passed through to the assembler.
        method: Assembler to use (``spades``, ``megahit`` or ``unicycler``).
        memory: Memory in GB for steps that take ``-m/--memory``; None keeps each step's default.
        minlen: Minimum read length after trimming.
        screen_accessions: GenBank accessions to screen out of the reads in ``filter``.
        screen_urls: URLs of sequences to screen out of the reads in ``filter``.
        mincontiglen: Minimum contig length kept by ``rmdup`` and ``sort``.
        workdir: Working directory for steps that take ``-w/--workdir``; None keeps each step's default.
        sourdb: ``taxpurge`` databases (names in ``SOURMASH_DATABASES``); None keeps its default.
        fcs_screen_run: Run ``fcs_screen`` (NCBI FCS-adaptor) after ``vecscreen``.
        container_engine: ``fcs_screen`` container engine; None keeps its default.
        mincovpct: Minimum percent of N50 coverage below which ``covpurge`` removes contigs.
        debug: Show debug output and keep temporary files in every step.
        quiet: Only show warnings and errors in every step.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ...); ignored.

    Raises:
        RuntimeError: If a step finishes without producing its expected output file.
    """
    fcs = _fcs_options(fcs_screen_run, container_engine)
    shared = _shared(cpus, memory, workdir, debug, quiet)
    filtered_1, filtered_2 = _trim_and_filter(read1, read2, basename, shared, minlen, screen_accessions, screen_urls)

    assembly = basename + f".{method}.fasta"
    _run_step(assemble, "assemble", assembly, shared, **_assemble_workdir(shared, method), read1=filtered_1, read2=filtered_2, out=assembly, method=method, tmpdir=tmpdir, assembler_args=assembler_args)

    _clean_and_finish(assembly, basename, shared, mincontiglen, _taxpurge_options(phylum, sourdb), {"read1": filtered_1, "read2": filtered_2, "mincovpct": mincovpct}, fcs=fcs, depth_options={"read1": filtered_1, "read2": filtered_2})


def run_long(
    longreads: str,
    basename: str,
    longread_type: str = "nano-hq",
    genome_size: str | None = None,
    phylum: list[str] | None = None,
    sourdb: list[str] | None = None,
    fcs_screen_run: bool = False,
    container_engine: str | None = None,
    cpus: int = 1,
    assembler_args: list[str] | None = None,
    memory: int | None = None,
    mincontiglen: int = 500,
    workdir: str | None = None,
    debug: bool = False,
    quiet: bool = False,
    **kwargs: Any,
) -> None:
    """Run the long-read (ONT/PacBio) pipeline, skipping steps whose output file already exists.

    ``taxpurge`` runs only when ``phylum`` is given; ``covpurge`` is not run, since it needs
    Illumina reads.

    Args:
        longreads: Long-read FASTQ.
        basename: Prefix for every output file (``{basename}.flye.fasta``, ``{basename}.final.fasta``, ...).
        longread_type: Flye read type (``nano-raw``, ``nano-hq``, ``pacbio-raw`` or ``pacbio-hifi``).
        genome_size: Estimated genome size for Flye (e.g. ``40m``), or None to let Flye estimate it.
        phylum: Phyla whose sourmash matches are kept by ``taxpurge``; None skips ``taxpurge``.
        sourdb: ``taxpurge`` databases (names in ``SOURMASH_DATABASES``); None keeps its default.
        fcs_screen_run: Run ``fcs_screen`` (NCBI FCS-adaptor) after ``vecscreen``.
        container_engine: ``fcs_screen`` container engine; None keeps its default.
        cpus: Threads for every step that takes ``-c/--cpus``.
        assembler_args: Extra arguments passed through to Flye.
        memory: Memory in GB for steps that take ``-m/--memory``; None keeps each step's default.
        mincontiglen: Minimum contig length kept by ``rmdup`` and ``sort``.
        workdir: Working directory for steps that take ``-w/--workdir``; None keeps each step's default.
        debug: Show debug output and keep temporary files in every step.
        quiet: Only show warnings and errors in every step.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ...); ignored.

    Raises:
        RuntimeError: If a step finishes without producing its expected output file.
    """
    fcs = _fcs_options(fcs_screen_run, container_engine)
    shared = _shared(cpus, memory, workdir, debug, quiet)
    assembly = _assemble_flye(longreads, basename, shared, longread_type, genome_size, assembler_args)
    depth_options = {"longreads": longreads, "longread_preset": LONGREAD_PRESETS[longread_type]}
    _clean_and_finish(assembly, basename, shared, mincontiglen, _taxpurge_options(phylum, sourdb), polish_steps=[{"method": "racon", "longreads": longreads}], fcs=fcs, depth_options=depth_options)


def run_hybrid(
    read1: str,
    longreads: str,
    basename: str,
    phylum: list[str] | None = None,
    read2: str | None = None,
    longread_type: str = "nano-hq",
    genome_size: str | None = None,
    method: str = "flye",
    cpus: int = 1,
    tmpdir: str | None = None,
    assembler_args: list[str] | None = None,
    memory: int | None = None,
    minlen: int = 75,
    screen_accessions: list[str] | None = None,
    screen_urls: list[str] | None = None,
    mincontiglen: int = 500,
    workdir: str | None = None,
    sourdb: list[str] | None = None,
    fcs_screen_run: bool = False,
    container_engine: str | None = None,
    mincovpct: int = 5,
    debug: bool = False,
    quiet: bool = False,
    **kwargs: Any,
) -> None:
    """Run the hybrid (Illumina + long-read) pipeline, skipping steps whose output file already exists.

    With ``method="flye"`` the long reads are assembled by Flye and, after the contigs are cleaned
    (vecscreen, taxpurge, covpurge, rmdup), polished by Racon, then with the filtered Illumina reads (Polypolish for paired reads, pypolca for single-end). With
    ``method="unicycler"`` Unicycler assembles both read types and polishes its own output, so the
    polish steps are skipped.

    Args:
        read1: Read 1 (forward, or single-end) raw Illumina FASTQ reads.
        longreads: Long-read FASTQ.
        basename: Prefix for every output file.
        phylum: Phyla whose sourmash matches are kept by ``taxpurge``; None skips ``taxpurge``.
        read2: Read 2 (reverse) raw Illumina FASTQ reads, or None for single-end data.
        longread_type: Flye read type (``nano-raw``, ``nano-hq``, ``pacbio-raw`` or ``pacbio-hifi``).
        genome_size: Estimated genome size for Flye (e.g. ``40m``), or None.
        method: ``flye`` or ``unicycler``.
        cpus: Threads for every step that takes ``-c/--cpus``.
        tmpdir: Assembler temporary directory.
        assembler_args: Extra arguments passed through to the assembler.
        memory: Memory in GB for steps that take ``-m/--memory``; None keeps each step's default.
        minlen: Minimum Illumina read length after trimming.
        screen_accessions: GenBank accessions to screen out of the Illumina reads in ``filter``.
        screen_urls: URLs of sequences to screen out of the Illumina reads in ``filter``.
        mincontiglen: Minimum contig length kept by ``rmdup`` and ``sort``.
        workdir: Working directory for steps that take ``-w/--workdir``; None keeps each step's default.
        sourdb: ``taxpurge`` databases (names in ``SOURMASH_DATABASES``); None keeps its default.
        fcs_screen_run: Run ``fcs_screen`` (NCBI FCS-adaptor) after ``vecscreen``.
        container_engine: ``fcs_screen`` container engine; None keeps its default.
        mincovpct: Minimum percent of N50 coverage below which ``covpurge`` removes contigs.
        debug: Show debug output and keep temporary files in every step.
        quiet: Only show warnings and errors in every step.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ...); ignored.

    Raises:
        ValueError: If ``method`` is not ``flye`` or ``unicycler``.
        RuntimeError: If a step finishes without producing its expected output file.
    """
    if method not in ("flye", "unicycler"):
        raise ValueError(f"pipeline_hybrid --method must be flye or unicycler, not {method}")
    fcs = _fcs_options(fcs_screen_run, container_engine)
    shared = _shared(cpus, memory, workdir, debug, quiet)
    filtered_1, filtered_2 = _trim_and_filter(read1, read2, basename, shared, minlen, screen_accessions, screen_urls)

    polish_steps: list[dict[str, Any]] = []
    if method == "unicycler":
        assembly = basename + ".unicycler.fasta"
        _run_step(assemble, "assemble", assembly, shared, **_assemble_workdir(shared, method), read1=filtered_1, read2=filtered_2, longreads=longreads, out=assembly, method=method, tmpdir=tmpdir, assembler_args=assembler_args)
    else:
        assembly = _assemble_flye(longreads, basename, shared, longread_type, genome_size, assembler_args)
        short_method = "polypolish" if filtered_2 else "pypolca"  # polypolish needs paired reads
        polish_steps = [{"method": "racon", "longreads": longreads}, {"method": short_method, "read1": filtered_1, "read2": filtered_2}]

    covpurge_options = {"read1": filtered_1, "read2": filtered_2, "mincovpct": mincovpct}
    _clean_and_finish(assembly, basename, shared, mincontiglen, _taxpurge_options(phylum, sourdb), covpurge_options, polish_steps, fcs=fcs, depth_options={"read1": filtered_1, "read2": filtered_2})


def _shared(cpus: int, memory: int | None, workdir: str | None, debug: bool, quiet: bool) -> dict[str, Any]:
    """Return the options passed to every step that has an option of the same name."""
    return {"cpus": cpus, "memory": memory, "workdir": workdir, "debug": debug, "quiet": quiet}


def _fcs_options(fcs_screen_run: bool, container_engine: str | None) -> dict[str, dict[str, Any] | None]:
    """Return the ``fcs_screen`` options (None when it is not run)."""
    screen = ({"container_engine": container_engine} if container_engine else {}) if fcs_screen_run else None
    return {"fcs_screen": screen}


def _taxpurge_options(phylum: list[str] | None, sourdb: list[str] | None) -> dict[str, Any] | None:
    """Return the ``taxpurge`` options, or None (skip ``taxpurge``) when no ``phylum`` is given.

    ``sourdb`` is passed only when given, so ``taxpurge`` otherwise keeps its own default database.
    """
    if not phylum:
        return None
    return {"phylum": phylum} | ({"sourdb": sourdb} if sourdb else {})


def _assemble_workdir(shared: dict[str, Any], method: str) -> dict[str, str]:
    """Return ``assemble``'s own ``workdir`` option: ``{workdir}/assemble_{method}`` when the pipeline has a workdir.

    The assemblers treat an existing output folder as a previous run of their own (SPAdes restarts
    from it, MEGAHIT refuses it), so they must not share the pipeline workdir that ``filter`` has
    already created. Without a pipeline workdir, ``assemble`` makes its own folder as usual.
    """
    return {"workdir": str(Path(shared["workdir"], f"assemble_{method}"))} if shared["workdir"] else {}


def _trim_and_filter(read1: str, read2: str | None, basename: str, shared: dict[str, Any], minlen: int, screen_accessions: list[str] | None, screen_urls: list[str] | None) -> tuple[str, str | None]:
    """Run ``trim`` then ``filter`` on the Illumina reads; return the filtered (read1, read2) files."""
    # trim writes {basename}_1P/_2P for paired reads, {basename}_1U for single-end reads
    trimmed_1, trimmed_2 = (basename + "_1P.fastq.gz", basename + "_2P.fastq.gz") if read2 else (basename + "_1U.fastq.gz", None)
    filtered_1, filtered_2 = basename + "_filtered_1.fastq.gz", (basename + "_filtered_2.fastq.gz" if read2 else None)
    _run_step(trim, "trim", trimmed_1, shared, read1=read1, read2=read2, basename=basename, minlen=minlen)
    _run_step(aaftf_filter, "filter", filtered_1, shared, read1=trimmed_1, read2=trimmed_2, basename=basename, screen_accessions=screen_accessions, screen_urls=screen_urls)
    return filtered_1, filtered_2


def _assemble_flye(longreads: str, basename: str, shared: dict[str, Any], longread_type: str, genome_size: str | None, assembler_args: list[str] | None) -> str:
    """Assemble ``longreads`` with Flye into ``{basename}.flye.fasta`` and return that file name."""
    assembly = basename + ".flye.fasta"
    _run_step(assemble, "assemble", assembly, shared, **_assemble_workdir(shared, "flye"), read1=None, longreads=longreads, out=assembly, method="flye", longread_type=longread_type, genome_size=genome_size, assembler_args=assembler_args)
    return assembly


def _clean_and_finish(
    assembly: str,
    basename: str,
    shared: dict[str, Any],
    mincontiglen: int,
    taxpurge_options: dict[str, Any] | None = None,
    covpurge_options: dict[str, Any] | None = None,
    polish_steps: list[dict[str, Any]] | tuple[()] = (),
    fcs: dict[str, dict[str, Any] | None] | None = None,
    depth_options: dict[str, Any] | None = None,
) -> None:
    """Run vecscreen → fcs_screen → taxpurge → covpurge → rmdup → polish → sort → depth → assess on ``assembly``.

    Polishing comes after the cleaning steps so contigs that are discarded are never polished.

    Args:
        assembly: Assembled FASTA.
        basename: Prefix for every output file.
        shared: Pipeline-wide options (see ``_step_kwargs``).
        mincontiglen: Minimum contig length kept by ``rmdup`` and ``sort``.
        taxpurge_options: ``taxpurge`` options (phylum, sourdb), or None to skip ``taxpurge``.
        covpurge_options: ``covpurge`` options (reads, mincovpct), or None to skip ``covpurge``.
        polish_steps: ``polish`` options for each polishing round, in order; each needs ``method``
            and writes ``{basename}.{method}.fasta``, which the next round (or ``sort``) reads.
        fcs: ``fcs_screen`` options (from ``_fcs_options``); skipped when None, or when ``fcs``
            is None.
        depth_options: ``depth`` read options; ``depth`` maps them to the sorted assembly and
            writes ``{basename}.final.fasta`` with likely mitochondrial contigs tagged (none is
            removed). When None, or when ``depth`` fails, the sorted assembly is used untagged.
    """
    vecscreen_file = basename + ".vecscreen.fasta"
    _run_step(vecscreen, "vecscreen", vecscreen_file, shared, infile=assembly, outfile=vecscreen_file)

    rmdup_input = vecscreen_file
    fcs = fcs or {}
    if fcs.get("fcs_screen") is not None:
        fcs_screen_file = basename + ".fcs_screen.fasta"
        _run_step(fcs_screen, "fcs_screen", fcs_screen_file, shared, infile=rmdup_input, outfile=fcs_screen_file, **fcs["fcs_screen"])
        rmdup_input = fcs_screen_file
    if taxpurge_options is not None:
        taxpurge_file = basename + ".taxpurge.fasta"
        _run_step(taxpurge, "taxpurge", taxpurge_file, shared, input=rmdup_input, outfile=taxpurge_file, **taxpurge_options)
        rmdup_input = taxpurge_file
    if covpurge_options is not None:
        covpurge_file = basename + ".covpurge.fasta"
        _run_step(covpurge, "covpurge", covpurge_file, shared, input=rmdup_input, outfile=covpurge_file, **covpurge_options)
        rmdup_input = covpurge_file

    rmdup_file = basename + ".rmdup.fasta"
    _run_step(rmdup, "rmdup", rmdup_file, shared, input=rmdup_input, out=rmdup_file, minlen=mincontiglen)

    polished = rmdup_file
    for options in polish_steps:
        polish_out = basename + f".{options['method']}.fasta"
        _run_step(polish, "polish", polish_out, shared, infile=polished, outfile=polish_out, **options)
        polished = polish_out

    sorted_file = basename + ".sorted.fasta"
    _run_step(aaftf_sort, "sort", sorted_file, shared, input=polished, out=sorted_file, minlen=mincontiglen)

    final_file = basename + ".final.fasta"
    _tag_mito(sorted_file, final_file, basename, shared, depth_options)

    assess.run(**_step_kwargs("assess", shared, input=final_file))


def _tag_mito(sorted_file: str, final_file: str, basename: str, shared: dict[str, Any], depth_options: dict[str, Any] | None) -> None:
    """Write ``final_file``: ``sorted_file`` with likely mitochondrial contigs tagged by ``depth``.

    ``depth`` also writes ``{basename}.coverage_stats.txt`` and ``{basename}.coverage_stats.mito_candidates.tsv``.
    Tagging is a convenience, so a failure only logs a warning and ``sorted_file`` is copied
    untagged; ``final_file`` is never missing because of it.
    """
    if check_file(final_file):
        logger.info(f"AAFTF depth output found: {final_file}")
        return
    if depth_options is not None:
        try:
            _run_step(depth, "depth", final_file, shared, input=sorted_file, out=basename + ".coverage_stats.txt", tagged=final_file, **depth_options)
            return
        except Exception as e:
            logger.warning(f"AAFTF depth failed, so mitochondrial contigs are not tagged: {e}")
    shutil.copyfile(sorted_file, final_file)


def _run_step(module: ModuleType, name: str, output: str, shared: dict[str, Any], **step_options: Any) -> None:
    """Run ``module.run()`` for step ``name`` unless ``output`` already exists, then check that it does.

    Args:
        module: The AAFTF subcommand module whose ``run()`` to call.
        name: Subcommand name, used to look up its CLI defaults and in log messages.
        output: File the step must produce; the step is skipped if it already exists.
        shared: Pipeline-wide options (see ``_step_kwargs``).
        **step_options: Step-specific keyword arguments that override the defaults.

    Raises:
        RuntimeError: If ``output`` is missing or empty after the step runs.
    """
    if check_file(output):
        logger.info(f"AAFTF {name} output found: {output}")
        return
    module.run(**_step_kwargs(name, shared, **step_options))
    if not check_file(output):
        raise RuntimeError(f"AAFTF {name} failed: {output} is missing or empty")


def _step_kwargs(name: str, shared: dict[str, Any], **step_options: Any) -> dict[str, Any]:
    """Return the run() keyword arguments for step ``name``: its CLI defaults, overridden by the pipeline's options.

    ``shared`` options are only passed to steps that have them, and a ``None``
    shared value (e.g. no pipeline ``--memory``) keeps the step's own default. ``pipe`` is always
    set to True so the steps don't log their "next command" hints (the pipeline runs the next step
    itself); it is a ``run()``-only parameter, not a CLI option.

    Args:
        name: Subcommand name whose CLI defaults to start from.
        shared: Pipeline-wide options (``cpus``, ``memory``, ``workdir``, ``debug``, ``quiet``).
        **step_options: Step-specific keyword arguments; these override everything else.

    Returns:
        Keyword arguments ready to pass to the step module's ``run()``.
    """
    kwargs = dict(_subcommand_defaults()[name])
    for key, value in shared.items():
        if key in kwargs and value is not None:
            kwargs[key] = value
    kwargs.update(step_options)
    kwargs["pipe"] = True
    return kwargs


@functools.cache
def _subcommand_defaults() -> dict[str, dict[str, Any]]:
    """Return ``{subcommand: {dest: default}}`` read from the AAFTF subcommand parsers."""
    from aaftf._menu import register_subcommands  # _menu imports this module, so import it when first needed

    subparsers = register_subcommands(ap.ArgumentParser())
    return {name: {action.dest: action.default for action in parser._actions if action.dest != "help"} for name, parser in subparsers.choices.items()}
