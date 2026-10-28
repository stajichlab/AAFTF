"""Utility scripts for parsing FASTA/FASTQ files and other shared helpers."""

import argparse as ap
import functools
import gzip
import logging
import os
import re
import shutil
import subprocess
import textwrap
import urllib.request
import uuid
from collections.abc import Callable, Iterable, Iterator, Sequence
from http.client import HTTPMessage
from itertools import islice
from pathlib import Path
from typing import IO, Any, BinaryIO, NamedTuple, TextIO, cast

import psutil
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from aaftf.resources import DATABASES

__all__ = [
    "COMPLEMENT",
    "CustomHelpFormatter",
    "PafHit",
    "SubcommandGroup",
    "bbtools_resources",
    "concat_files",
    "dedupe_fasta",
    "open_maybe_gz",
    "estimate_read_length",
    "check_file",
    "available_cpus",
    "get_ram",
    "fasta_stats",
    "filter_fasta",
    "write_fasta",
    "softwrap",
    "count_fastq",
    "align_to_sorted_bam",
    "samtools_sort_cmd",
    "print_cmd",
    "bam_read_count",
    "paf_hits",
    "execute",
    "calc_nx",
    "run_cmd",
    "require_tools",
    "next_step_name",
    "basename_from_reads",
    "require_databases",
    "find_db_file",
    "find_db_index",
    "db_dirs",
    "home_db_cache",
    "db_file",
    "db_write_dir",
    "download_file",
    "warn_if_home_cache",
    "open_url",
    "safe_remove",
    "setup_logging",
    "finish_logging",
    "make_workdir",
    "cleanup_workdir",
    "close_step_log",
]


# IUPAC complement (U pairs with A); case is preserved so soft-masked bases stay lowercase.
COMPLEMENT = str.maketrans("ACGTURYKMSWBDHVNXacgturykmswbdhvnx", "TGCAAYRMKSWVBDHNXtgcaayrmkswvbdhnx")

logger = logging.getLogger(__name__)

_home_cache_warned = False
# per-step log file state (see setup_logging / make_workdir / cleanup_workdir / finish_logging)
_log_buffer: "_BufferHandler | None" = None
_step_log: logging.FileHandler | None = None
_step_log_used = False
_last_step_log: Path | None = None


class CustomHelpFormatter(ap.HelpFormatter):
    """Custom help formatter for argparse to enhance text wrapping and default value display.

    This formatter adjusts the text wrapping for better readability and ensures that the default value of an argument is displayed
    in the help message if not already present.
    """

    def _fill_text(self, text: str, width: int, indent: str) -> str:
        """Re-wrap a description/epilog paragraph by paragraph, keeping blank-line breaks.

        Args:
            text: Description or epilog text; paragraphs are separated by blank lines.
            width: Maximum line width.
            indent: Indentation argparse asks for; ignored.

        Returns:
            The re-wrapped text.
        """
        paragraphs = [self._whitespace_matcher.sub(" ", paragraph.strip()) for paragraph in text.split("\n\n") if paragraph.strip()]
        return "\n\n".join([textwrap.fill(paragraph, width) for paragraph in paragraphs])

    def _split_lines(self, text: str, width: int) -> list[str]:
        """Wrap argument help, keeping the help string's own line breaks.

        Args:
            text: Argument help text.
            width: Maximum line width.

        Returns:
            The wrapped lines.
        """
        lines = [self._whitespace_matcher.sub(" ", line.strip()) for line in text.split("\n") if line.strip()]
        return [wrapped_line for line in lines for wrapped_line in textwrap.wrap(line, width)]

    def _get_help_string(self, action: ap.Action) -> str | None:
        """Append ``(default: ...)`` to an argument's help unless it already mentions its default.

        Nothing is appended when the default is suppressed, None or False, or when the argument
        has no help text.

        Args:
            action: The argument whose help is being formatted.

        Returns:
            The help string (a %-format template argparse fills in).
        """
        help = action.help
        if help is None:
            return None
        pattern = r"\(default: .+\)"
        if re.search(pattern, help) is None:
            if action.default is not ap.SUPPRESS and action.default is not None and action.default is not False:
                defaulting_nargs = [ap.OPTIONAL, ap.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help

    def _format_action(self, action: ap.Action) -> str:
        """List a SubcommandGroup's subcommands directly under its group title, without a header line.

        Args:
            action: The action to format; any non-``SubcommandGroup`` is formatted as usual.

        Returns:
            The formatted help text for the action.
        """
        if isinstance(action, SubcommandGroup):
            return "".join(self._format_action(sub) for sub in action.subcommands)
        return super()._format_action(action)


class SubcommandGroup(ap.Action):
    """Help-only entry listing some of a parser's subcommands under their own argument-group title.

    argparse prints every subcommand in one block; add one of these to each argument group
    (``group._group_actions.append(SubcommandGroup(...))``) and hide the real subparsers action
    (``help=argparse.SUPPRESS``) to show them in sections. It is never parsed.
    """

    def __init__(self, subcommands: list[ap.Action]) -> None:
        """Store the subcommand help entries to list.

        Args:
            subcommands: Help entries taken from the subparsers action's ``_choices_actions``.
        """
        super().__init__(option_strings=[], dest=ap.SUPPRESS, nargs=0, metavar="")
        self.subcommands = subcommands

    def _get_subactions(self) -> list[ap.Action]:
        """Let the help formatter size its columns to fit the listed subcommands.

        Returns:
            The listed subcommand help entries.
        """
        return self.subcommands

    def __call__(self, parser: ap.ArgumentParser, namespace: ap.Namespace, values: Any, option_string: str | None = None) -> None:
        """Do nothing; never called, because this action is not registered with the parser.

        Args:
            parser: Unused.
            namespace: Unused.
            values: Unused.
            option_string: Unused.
        """


class PafHit(NamedTuple):
    """One alignment from minimap2's PAF output (its 12 standard columns; coordinates are 0-based).

    Attributes:
        query: Query sequence name.
        query_len: Query sequence length.
        query_start: Query start (0-based).
        query_end: Query end (exclusive).
        strand: ``+`` or ``-``.
        target: Target sequence name.
        target_len: Target sequence length.
        target_start: Target start (0-based).
        target_end: Target end (exclusive).
        matches: Number of matching bases.
        aln_len: Alignment block length, including gaps.
        mapq: Mapping quality (0-255).
    """

    query: str
    query_len: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_len: int
    target_start: int
    target_end: int
    matches: int
    aln_len: int
    mapq: int


def concat_files(paths: Iterable[str | Path], dest: str | Path) -> None:
    """Concatenate files (gzipped ones are decompressed) into ``dest``.

    Args:
        paths: Files to concatenate, in order.
        dest: Output file (overwritten).
    """
    with open(dest, "wb") as out:
        for path in paths:
            with open_maybe_gz(path, "rb") as fh:
                shutil.copyfileobj(fh, out)


def open_maybe_gz(path: str | Path, mode: str = "rt") -> IO[Any]:
    """Open ``path`` with gzip if its name ends in ``.gz``, else as a plain file.

    Args:
        path: File to open.
        mode: Open mode, e.g. ``"rt"`` or ``"rb"``.

    Returns:
        The open file object.
    """
    return cast(IO[Any], (gzip.open if str(path).endswith(".gz") else open)(path, mode))


def estimate_read_length(input: str | Path) -> int:
    """Guess the read length in a FASTQ file, rounded to the nearest 10bp.

    Uses the longest of the first 500 reads.

    Args:
        input: FASTQ file, optionally gzipped.

    Returns:
        The estimated read length.
    """
    with open_maybe_gz(input) as infile:
        records = islice(FastqGeneralIterator(infile), 500)
        max_len = max(len(seq) for _, seq, _ in records)
    return round(max_len, -1)


def check_file(input: str | Path) -> bool:
    """Check that ``input`` is an existing, non-empty regular file.

    Args:
        input: File path to check.

    Returns:
        True if the file exists and is at least 1 byte long.
    """

    def _get_size(filename: str | Path) -> int:
        """Return the size of ``filename`` in bytes.

        Args:
            filename: File path.

        Returns:
            File size in bytes.
        """
        st = os.stat(filename)
        return st.st_size

    if Path(input).is_file():
        return _get_size(input) >= 1
    return False


def available_cpus() -> int:
    """Return how many CPUs this process may use.

    Checks, in order: the SLURM allocation (``SLURM_CPUS_PER_TASK``), the CPU
    affinity set (respects ``taskset``, SLURM CPU binding and cgroup cpusets),
    then the machine's CPU count where affinity is unavailable (macOS, Windows).
    CPU time quotas such as ``docker --cpus`` are not detected.

    Returns:
        Number of usable CPUs (at least 1).
    """
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK", "")
    if slurm_cpus.isdigit() and int(slurm_cpus) > 0:
        return int(slurm_cpus)
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def _read_int(path: str) -> int | None:
    """Return the integer in a one-value sysfs file, or None if missing or not a number (e.g. ``max``)."""
    try:
        with open(path) as fh:
            return int(fh.read().strip())
    except (OSError, ValueError):
        return None


def _inactive_file(stat_path: str, key: str) -> int:
    """Return the reclaimable page-cache bytes (``key``) from a cgroup ``memory.stat`` file, or 0."""
    try:
        with open(stat_path) as fh:
            for line in fh:
                name, _, value = line.partition(" ")
                if name == key:
                    return int(value)
    except (OSError, ValueError):
        pass
    return 0


def _cgroup_memory_available(proc_cgroup: str = "/proc/self/cgroup", sysfs: str = "/sys/fs/cgroup") -> int | None:
    """Return the bytes still free under the tightest cgroup memory limit, or None if there is none.

    SLURM (``ConstrainRAMSpace``), systemd user slices on login nodes and
    containers limit memory with cgroups, which ``/proc/meminfo`` does not
    show. The limit is usually set on a parent cgroup (e.g. ``job_<id>``, not
    the ``task_0`` the process sits in), so every cgroup from the process's own
    up to the root is checked. Free = limit - (usage - reclaimable page cache).

    Args:
        proc_cgroup: The process's cgroup membership file.
        sysfs: Mount point of the cgroup filesystem.

    Returns:
        Free bytes under the tightest limit (at least 0), or None when unlimited or unknown.
    """
    try:
        with open(proc_cgroup) as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None
    best = None
    for line in lines:
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        _, controllers, path = parts
        if controllers == "":  # cgroup v2
            root = sysfs
            limit_f, usage_f, cache_key = "memory.max", "memory.current", "inactive_file"
        elif "memory" in controllers.split(","):  # cgroup v1
            root = os.path.join(sysfs, "memory")
            limit_f, usage_f, cache_key = "memory.limit_in_bytes", "memory.usage_in_bytes", "total_inactive_file"
        else:
            continue
        d = os.path.normpath(root + "/" + path.lstrip("/"))
        while True:
            limit = _read_int(os.path.join(d, limit_f))
            if limit is not None and limit < 1 << 60:  # v1 reports "unlimited" as a huge number
                usage = _read_int(os.path.join(d, usage_f)) or 0
                usage -= _inactive_file(os.path.join(d, "memory.stat"), cache_key)
                free = max(limit - max(usage, 0), 0)
                best = free if best is None else min(best, free)
            if d == root or not d.startswith(root):
                break
            d = os.path.dirname(d)
    return best


def get_ram(max_lim: float = 0) -> float:
    """Get the RAM available to this process, in GB, kept safely under the true value.

    The smallest of: the machine's available memory (``psutil``), the free
    memory under any cgroup limit (SLURM job, login-node user slice,
    container), the SLURM allocation (``SLURM_MEM_PER_NODE``) and the process
    address-space limit (``ulimit -v``). Off Linux only the first applies.

    Rounded down to the nearest 10 once available RAM exceeds 10 GB;
    otherwise rounded down to 1 decimal with a small safety margin
    subtracted, clamped at 0. This never exceeds the true available RAM.

    Args:
        max_lim: Optional upper bound (GB) to cap the result at; 0 means no cap.

    Returns:
        Available RAM in GB.
    """
    avail = psutil.virtual_memory().available
    cgroup_free = _cgroup_memory_available()
    if cgroup_free is not None:
        avail = min(avail, cgroup_free)
    slurm_mem = os.environ.get("SLURM_MEM_PER_NODE", "")
    if slurm_mem.isdigit() and int(slurm_mem) > 0:
        avail = min(avail, int(slurm_mem) * 1024**2)  # MB
    rlimit_as = getattr(psutil, "RLIMIT_AS", None)
    if rlimit_as is not None:
        try:
            soft, _ = psutil.Process().rlimit(rlimit_as)
            if soft != psutil.RLIM_INFINITY and soft > 0:
                avail = min(avail, soft)
        except (OSError, psutil.Error):
            pass
    avail_gb = avail / (1024.0**3)
    safe_gb = avail_gb // 10 * 10 if avail_gb > 10 else max(round(avail_gb, 1) - 0.1, 0)
    return min(safe_gb, max_lim) if max_lim else safe_gb


def fasta_stats(input: str | Path) -> tuple[int, int]:
    """Return (number of records, total sequence length) of a FASTA file.

    Args:
        input: Uncompressed FASTA file.

    Returns:
        Tuple (number of records, total sequence length).
    """
    count = length = 0
    with open(input) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
            else:
                length += len(line.rstrip())
    return count, length


def filter_fasta(fasta_in: str | Path, fasta_out: str | Path, keep: Callable[[str], Any], wrap: int = 60) -> tuple[int, int]:
    """Write records whose ID passes ``keep(id)`` to fasta_out, wrapped at ``wrap`` columns.

    The ID is the first whitespace-delimited word of the header (as Biopython's
    ``record.id``); the full header line is kept. The default 60-column wrap
    matches Biopython's ``SeqIO.write``.

    Args:
        fasta_in: Input FASTA file.
        fasta_out: Output FASTA file (overwritten).
        keep: Predicate called with each record ID; truthy keeps the record.
        wrap: Sequence line width.

    Returns:
        Tuple (number of records written, total length written).
    """
    count = length = 0
    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        for header, seq in SimpleFastaParser(fin):
            if keep(header.split(None, 1)[0]):
                write_fasta(fout, header, seq, wrap)
                count += 1
                length += len(seq)
    return count, length


def write_fasta(fh: TextIO, header: str, seq: str, wrap: int = 60) -> None:
    """Write one FASTA record to an open file handle, wrapped at ``wrap`` columns.

    Args:
        fh: Text file handle open for writing.
        header: Header line without the leading ``>``.
        seq: Sequence; if empty, only the header is written.
        wrap: Sequence line width.
    """
    fh.write(f">{header}\n")
    if seq:
        fh.write(f"{softwrap(seq, wrap)}\n")


def softwrap(string: str, every: int = 60) -> str:
    """Softwrap lines in a textstring (60 columns by default, as Biopython's SeqIO.write).

    Args:
        string: Text to wrap.
        every: Line width.

    Returns:
        The text split into newline-joined chunks of ``every`` characters.
    """
    return "\n".join(string[i : i + every] for i in range(0, len(string), every))


def count_fastq(input: str) -> int:
    """Count the number of records in a FASTQ file (gzip or regular).

    Gzipped input is decompressed with pigz (or gzip) when available, which is
    much faster than Python's gzip module.

    Args:
        input: FASTQ file path; a ``.gz`` suffix marks it as gzipped.

    Returns:
        Number of records (line count divided by 4).

    Raises:
        OSError: If the file is missing or not valid gzip.
        subprocess.CalledProcessError: If pigz/gzip fails to decompress it.
    """
    decompressor = input.endswith(".gz") and (shutil.which("pigz") or shutil.which("gzip"))
    if decompressor:
        with subprocess.Popen([decompressor, "-dc", input], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as proc:
            assert proc.stdout is not None  # stdout=PIPE always sets it
            lines = _count_lines(proc.stdout)
        if proc.returncode:
            raise subprocess.CalledProcessError(proc.returncode, [decompressor, "-dc", input])
    else:
        with open_maybe_gz(input, "rb") as fh:
            lines = _count_lines(fh)
    return lines // 4


def align_to_sorted_bam(align_cmd: Sequence[str], bam_out: str | Path, threads: int = 1, cwd: str | Path | None = None, debug: bool = False) -> None:
    """Pipe an aligner's SAM output straight into ``samtools sort``.

    The BAM is indexed (``<bam_out>.bai``) as it is written. If either command
    fails, any partial BAM/index is removed so a rerun does not mistake it for
    a finished alignment.

    Args:
        align_cmd: Aligner command list that writes SAM to stdout.
        bam_out: Destination sorted BAM, relative to the current directory
            (not ``cwd``).
        threads: samtools sort threads.
        cwd: Working directory for both commands.
        debug: Show both commands' stderr (hidden otherwise).

    Raises:
        RuntimeError: If the aligner or ``samtools sort`` exits non-zero.
    """
    bam_path = Path(bam_out).resolve()
    sort_cmd = samtools_sort_cmd("-", str(bam_path), threads, write_index=True)
    stderr = None if debug else subprocess.DEVNULL
    print_cmd(align_cmd)
    print_cmd(sort_cmd)
    p1 = subprocess.Popen(align_cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=stderr)
    p2 = subprocess.Popen(sort_cmd, cwd=cwd, stdin=p1.stdout, stderr=stderr)
    assert p1.stdout is not None  # stdout=PIPE always sets it
    p1.stdout.close()
    p2.communicate()
    if p1.wait() != 0 or p2.returncode != 0:
        bam_path.unlink(missing_ok=True)
        Path(f"{bam_path}.bai").unlink(missing_ok=True)
        raise RuntimeError(f"{align_cmd[0]} | samtools sort failed for {bam_out}")


def samtools_sort_cmd(input_file: str, output_bam: str, threads: int = 1, memory_per_thread: str | None = None, tmp_prefix: str | None = None, write_index: bool = False) -> list[str]:
    """Return a ``samtools sort`` command list (samtools >= 1.13; pixi pins >= 1.24).

    Args:
        input_file: Path to input SAM/BAM, or ``'-'`` to read from stdin.
        output_bam: Destination sorted BAM path.
        threads: Number of sort threads.
        memory_per_thread: Memory per thread string, e.g. ``'1G'``.
        tmp_prefix: Prefix for temporary sort files (passed as -T).
        write_index: Also write ``<output_bam>.bai`` while sorting.

    Returns:
        The command as a list of arguments.
    """
    cmd = ["samtools", "sort", "-@", str(threads)]
    if memory_per_thread:
        cmd += ["-m", str(memory_per_thread)]
    if write_index:
        cmd += ["--write-index", "-o", f"{output_bam}##idx##{output_bam}.bai"]
    else:
        cmd += ["-o", output_bam]
    if tmp_prefix:
        cmd += ["-T", tmp_prefix]
    cmd.append(input_file)
    return cmd


def print_cmd(cmd: Sequence[str]) -> None:
    """Print a command, wrapped at 80 columns with a coloured ``CMD:`` prefix.

    Args:
        cmd: Command list.
    """
    stringcmd = "{:}".format(" ".join(cmd))
    prefix = "\033[96mCMD:\033[00m "
    wrapper = textwrap.TextWrapper(initial_indent=prefix, width=80, subsequent_indent=" " * 8, break_long_words=False)
    print(wrapper.fill(stringcmd))


def bam_read_count(bamfile: str) -> tuple[int, int]:
    """Return (mapped, unmapped) counts of reads in a BAM file.

    Uses ``samtools flagstat`` primary counts, so secondary and supplementary
    alignments are not counted as extra reads.

    Args:
        bamfile: BAM file path.

    Returns:
        Tuple (mapped reads, unmapped reads).
    """
    stats = {}
    for line in execute(["samtools", "flagstat", "-O", "tsv", bamfile], quiet=True):
        passed, _failed, label = line.rstrip("\n").split("\t")
        stats[label] = passed
    primary, mapped = int(stats["primary"]), int(stats["primary mapped"])
    return mapped, primary - mapped


# from https://stackoverflow.com/questions/4417546/
# constantly-print-subprocess-output-while-process-is-running
def paf_hits(cmd: Sequence[str], **execute_args: Any) -> Iterator[PafHit]:
    """Run a PAF-producing command (e.g. ``minimap2 -x ...``) and yield each alignment as a ``PafHit``.

    Args:
        cmd: Command list whose stdout is PAF.
        **execute_args: Passed to ``execute`` (``cwd``, ``debug``, ``quiet``).

    Yields:
        One ``PafHit`` per PAF line with at least 12 columns.
    """
    for line in execute(cmd, **execute_args):
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 12:
            continue
        yield PafHit(cols[0], int(cols[1]), int(cols[2]), int(cols[3]), cols[4], cols[5], *map(int, cols[6:12]))


def execute(cmd: Sequence[str], cwd: str | Path | None = None, debug: bool = False, quiet: bool = False) -> Iterator[str]:
    """Run a command and yield its stdout line by line.

    Args:
        cmd: Command list.
        cwd: Working directory.
        debug: Show the command's stderr (hidden otherwise).
        quiet: Don't print the command (e.g. when it runs once per contig).

    Yields:
        Each stdout line, including its trailing newline.

    Raises:
        subprocess.CalledProcessError: If the command fails and all its output
            was read. If the caller stops reading early, the command is
            stopped and its exit status is ignored.
    """
    if not quiet:
        print_cmd(cmd)
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, text=True, stderr=None if debug else subprocess.DEVNULL)
    assert proc.stdout is not None  # stdout=PIPE always sets it
    finished = False
    try:
        yield from proc.stdout
        finished = True
    finally:
        proc.stdout.close()
        if not finished:
            proc.kill()
        return_code = proc.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def calc_nx(lengths: Sequence[int], fraction: float = 0.5) -> tuple[int, int]:
    """Return (NX, LX) for a set of contig lengths.

    NX is the length of the contig at which the longest contigs first cover
    ``fraction`` of the total length; LX is how many contigs that takes.
    The input is not modified.

    Args:
        lengths: Contig lengths.
        fraction: Fraction of the total length to cover, e.g. 0.5 for N50.

    Returns:
        Tuple (NX, LX); (0, 0) for an empty input.
    """
    target = sum(lengths) * fraction
    cumulative = 0
    for count, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative >= target:
            return length, count
    return 0, 0


def run_cmd(cmd: Sequence[str], debug: bool = False, cwd: str | Path | None = None, stdout: IO[Any] | int | None = None, env: dict[str, str] | None = None, quiet_stdout: bool = False) -> subprocess.CompletedProcess[bytes]:
    """Print a command, then run it with stderr hidden unless ``debug``.

    Args:
        cmd: Command list.
        debug: Show the command's stderr (and stdout, if ``quiet_stdout``).
        cwd: Working directory.
        stdout: stdout destination (e.g. an open file); inherited when None.
        env: Environment for the command; inherited when None.
        quiet_stdout: Also hide stdout unless ``debug``.

    Returns:
        The finished process; the exit status is not checked.
    """
    print_cmd(cmd)
    if quiet_stdout and not debug:
        stdout = subprocess.DEVNULL
    return subprocess.run(cmd, cwd=cwd, stdout=stdout, stderr=None if debug else subprocess.DEVNULL, env=env)


def bbtools_resources(names: Iterable[str]) -> list[str]:
    """Return the paths of files ``names`` in the BBTools ``resources/`` folder.

    The folder sits next to the real ``bbduk.sh`` script (after resolving the symlink conda/pixi
    put on ``PATH``).

    Args:
        names: File names in the resources folder, e.g. ``["adapters.fa"]``.

    Returns:
        Paths of the files, in the order of ``names``.

    Raises:
        FileNotFoundError: If ``bbduk.sh`` is not on ``PATH`` or the files are not in its resources folder.
    """
    bbduk = shutil.which("bbduk.sh")
    if not bbduk:
        raise FileNotFoundError("bbduk.sh not found on PATH; BBTools (bbmap) provides the bundled adapter and lambda phage sequences")
    resources = Path(bbduk).resolve().parent / "resources"
    paths = [resources / name for name in names]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"BBTools resource file(s) not found: {', '.join(missing)}")
    return [str(path) for path in paths]


def bbduk_reads(bbduk_cmd: Sequence[str], read1: str, read2: str | None, out1: str, out2: str | None, tmp_prefix: str | Path, debug: bool = False) -> None:
    """Run ``bbduk_cmd`` (``bbduk.sh`` and its options, no in/out) on the reads.

    Paired reads are interleaved with ``reformat.sh``, run through BBDuk's single-file reader
    and split back into ``out1``/``out2``: BBDuk's paired mode (``in=``/``in2=``, still in
    BBTools 39.91) crashes ("List size mismatch") or silently truncates its output when the
    two reads of a pair differ in length. The interleaved files are removed afterwards.

    Args:
        bbduk_cmd: The ``bbduk.sh`` command without input/output arguments.
        read1: Read 1 (forward, or single-end) FASTQ.
        read2: Read 2 (reverse) FASTQ, or None for single-end.
        out1: Output read 1 (or single-end) FASTQ.
        out2: Output read 2 FASTQ (used with ``read2``).
        tmp_prefix: Path prefix of the temporary interleaved files.
        debug: Show external command output when True.
    """
    if not read2:
        run_cmd([*bbduk_cmd, f"in={read1}", f"out={out1}"], debug)
        return
    interleaved_in = f"{tmp_prefix}_ivl.fq.gz"
    interleaved_out = f"{tmp_prefix}_ivl.bbduk.fq.gz"
    try:
        run_cmd(["reformat.sh", f"in1={read1}", f"in2={read2}", f"out={interleaved_in}"], debug)
        run_cmd([*bbduk_cmd, f"in={interleaved_in}", "interleaved=true", f"out={interleaved_out}"], debug)
        run_cmd(["reformat.sh", f"in={interleaved_out}", f"out1={out1}", f"out2={out2}"], debug)
    finally:
        safe_remove(interleaved_in)
        safe_remove(interleaved_out)


# dedupe.sh seed length: the shortest contamdb sequences (adaptors/primers) are 10 bp
DEDUPE_K = 10


def dedupe_fasta(src: str | Path, dest: str | Path, debug: bool = False) -> int:
    """Write ``src`` to ``dest`` without duplicate or contained sequences, with BBTools ``dedupe.sh``.

    Exact duplicates and sequences fully contained in another, in either orientation, are removed
    at 100% identity (an ``N`` matches only ``N``); partial overlaps are not merged. The seed
    length is ``DEDUPE_K`` rather than dedupe.sh's default 31, since sequences shorter than it are
    never found contained and adaptor/primer entries can be as short as 10 bp.

    Args:
        src: Input FASTA.
        dest: Output FASTA (overwritten); its name must end in a FASTA extension (e.g. ``.fa``).
        debug: Show ``dedupe.sh``'s log.

    Returns:
        Number of sequences written.

    Raises:
        FileNotFoundError: If ``dedupe.sh`` is not on PATH.
        RuntimeError: If ``dedupe.sh`` fails or writes no output.
    """
    require_tools(["dedupe.sh"], hint="BBTools (bbmap) provides dedupe.sh")
    cmd = ["dedupe.sh", f"in={src}", f"out={dest}", f"k={DEDUPE_K}", "absorbmatch=t", "absorbcontainment=t", "absorbrc=t", "exact=t", "findoverlap=f", "ow=t"]
    result = run_cmd(cmd, debug)
    if result.returncode != 0 or not Path(dest).is_file():
        raise RuntimeError(f"dedupe.sh failed (exit status {result.returncode}) on {src}")
    with open(dest) as fh:
        return sum(line.startswith(">") for line in fh)


def require_tools(tools: Iterable[str], hint: str | None = None) -> None:
    """Raise FileNotFoundError naming any of ``tools`` that are not on PATH.

    Args:
        tools: Executable names to look up.
        hint: Advice appended to the error; a generic install hint when None.

    Raises:
        FileNotFoundError: If any tool is missing from PATH.
    """
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        hint = hint or "Install them (e.g. `conda install -c bioconda <tool>`) and make sure the correct environment is activated."
        raise FileNotFoundError(f"required tool(s) not found on PATH: {', '.join(missing)}\n{hint}")


def next_step_name(outfile: str, suffix: str) -> str:
    """Suggest the next step's output name: ``outfile`` up to its first '_' (else first '.') plus ``suffix``.

    Args:
        outfile: This step's output filename.
        suffix: Suffix for the next step, e.g. ``".taxpurge.fasta"``.

    Returns:
        The suggested filename.
    """
    for sep in ("_", "."):
        if sep in outfile:
            return outfile.split(sep)[0] + suffix
    return outfile + suffix


def basename_from_reads(reads: str | Path) -> str:
    """Derive a sample basename from a reads filename: the name up to its first '_' (else first '.').

    Args:
        reads: Reads file path; only its final component is used.

    Returns:
        The sample basename.
    """
    name = Path(reads).name
    for sep in ("_", "."):
        if sep in name:
            return name.split(sep)[0]
    return name


def require_databases(names: Iterable[str], hint: str | None = None) -> list[str]:
    """Return the stored paths of catalog databases ``names`` (keys of ``resources.DATABASES``).

    Subcommands never download these themselves: if any is missing from every
    ``db_dirs()`` folder, raise an error naming the ``AAFTF database`` command
    that fetches them.

    Args:
        names: Database abbreviations, e.g. ``["phix", "euk_adaptors"]``.
        hint: Extra alternative to suggest, e.g. "or pass --sourdb PATH".

    Returns:
        Paths of the database files, in the order of ``names``.

    Raises:
        FileNotFoundError: If any database is not found.
    """
    paths, missing = [], []
    for name in names:
        path = find_db_file(DATABASES[name]["filename"])
        if path:
            paths.append(path)
        else:
            missing.append(name)
    if missing:
        searched = ", ".join(str(folder) for folder in db_dirs())
        suggestion = f"Download them first: AAFTF database {' '.join(missing)}" + (f" ({hint})" if hint else "")
        raise FileNotFoundError(f"missing database(s): {', '.join(missing)} (searched {searched})\n{suggestion}")
    return paths


def find_db_file(name: str) -> str | None:
    """Return the first existing copy of database file ``name`` in ``db_dirs()``, or None.

    A copy that differs only by a ``.gz`` suffix also matches (e.g. a user who
    gunzipped ``contam_in_euks.fa.gz``, or gzipped ``UniVec``); within each
    folder the exact name is preferred.

    Args:
        name: Database filename.

    Returns:
        Path of the first copy found, or None.
    """
    alt = name[: -len(".gz")] if name.endswith(".gz") else f"{name}.gz"
    for folder in db_dirs():
        for candidate in (name, alt):
            if (folder / candidate).is_file():
                return str(folder / candidate)
    return None


def find_db_index(name: str) -> str | None:
    """Return the first RocksDB index folder ``name`` in ``db_dirs()``, or None.

    Only a folder holding RocksDB's ``CURRENT`` file counts, so an empty or partly copied folder
    is not used.

    Args:
        name: Index folder name, e.g. ``gtdb-rs226-reps.k31.rocksdb``.

    Returns:
        Path of the first index found, or None.
    """
    for folder in db_dirs():
        if (folder / name / "CURRENT").is_file():
            return str(folder / name)
    return None


def db_dirs() -> list[Path]:
    """Return the database folders to search, in order.

    ``$AAFTF_DB`` may list several folders separated like ``$PATH`` (``:`` on
    Linux/macOS), e.g. a shared read-only folder followed by a personal one.
    When it is unset, the home cache (``home_db_cache()``) is used.

    Returns:
        Resolved folder paths.
    """
    folders = [Path(p).expanduser().resolve() for p in os.environ.get("AAFTF_DB", "").split(os.pathsep) if p]
    return folders or [home_db_cache()]


def home_db_cache() -> Path:
    """Return the default database folder: ``$XDG_CACHE_HOME/aaftf``, else ``~/.cache/aaftf``.

    Returns:
        The resolved folder path (not created).
    """
    return (Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "aaftf").expanduser().resolve()


def db_file(name: str, force: bool = False) -> str:
    """Return the path to use for database file ``name``.

    An existing copy in any ``db_dirs()`` folder is reused; otherwise (or when
    ``force`` asks for a fresh download) it is placed in ``db_write_dir()``.

    Args:
        name: Database filename.
        force: Ignore existing copies.

    Returns:
        The file path.
    """
    existing = None if force else find_db_file(name)
    return existing or str(db_write_dir() / name)


def db_write_dir() -> Path:
    """Return the folder new databases are downloaded to (created if needed).

    This is the first ``db_dirs()`` folder that can be written to, falling back
    to the home cache if none can.

    Returns:
        The writable folder.
    """
    for folder in db_dirs():
        try:
            folder.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        if os.access(folder, os.W_OK):
            return folder
    fallback = home_db_cache()
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def download_file(url: str, dest: str, force: bool = False) -> str:
    """Download ``url`` to ``dest`` unless it already exists.

    Writes to a temporary file and renames it on success, so an interrupted
    download never leaves a partial file that a later run would reuse.

    Args:
        url: Remote URL to download.
        dest: Local file path to write.
        force: Re-download even if ``dest`` exists.

    Returns:
        ``dest``.

    Raises:
        Exception: Any download or write error is logged and re-raised after the
            temporary file is removed.
    """
    if Path(dest).exists() and not force:
        logger.info(f"Already present: {dest}")
        return dest

    if Path(dest).expanduser().resolve().parent == home_db_cache():
        warn_if_home_cache()
    logger.info(f"Downloading {Path(dest).name} ...")
    Path(dest).parent.mkdir(parents=True, exist_ok=True)

    tmp = f"{dest}.tmp"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "AAFTF/1.0"})
        with open_url(req, timeout=300) as response:
            final_url = response.geturl()
            if final_url != url:
                logger.info(f"Redirected to {final_url}")
            with open(tmp, "wb") as outfh:
                shutil.copyfileobj(response, outfh)
        os.replace(tmp, dest)
    except Exception as e:
        logger.error(f"downloading {url}: {e}")
        safe_remove(tmp)
        raise

    logger.info(f"Saved {dest}")
    return dest


def warn_if_home_cache() -> None:
    """Warn once if databases are going to the home cache because ``$AAFTF_DB`` is unset or unwritable."""
    global _home_cache_warned
    home = home_db_cache()
    user_chose_home = bool(os.environ.get("AAFTF_DB")) and home in db_dirs()
    if _home_cache_warned or user_chose_home or db_write_dir() != home:
        return
    _home_cache_warned = True
    reason = "no folder in AAFTF_DB is writable" if os.environ.get("AAFTF_DB") else "AAFTF_DB is not set"
    logger.warning(f"{reason}, so databases are saved to {home}. Some are large (the sourmash databases and the FCS container image are several GB each) and may exceed a home-directory quota. To store them elsewhere: export AAFTF_DB=/path/with/space")


def open_url(request: str | urllib.request.Request, timeout: float) -> Any:
    """Open a URL or ``urllib.request.Request``, following HTTP 308 redirects (which Python < 3.11 does not).

    Args:
        request: URL string or prepared request.
        timeout: Socket timeout in seconds.

    Returns:
        The response object (usable as a context manager).
    """
    return _url_opener().open(request, timeout=timeout)


def safe_remove(path: str | Path) -> None:
    """Remove a file, symlink or directory tree; do nothing if it does not exist.

    A symlink is removed itself, never the directory it points to.

    Args:
        path: Path to remove.
    """
    path = Path(path)
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def setup_logging(debug: bool = False, quiet: bool = False) -> None:
    """Send AAFTF log messages to stderr, and hold them until a step log file is opened.

    The terminal shows INFO and above (only warnings and errors with ``quiet``; DEBUG too with
    ``debug``). Every message at INFO and above (DEBUG too with ``debug``) is also kept in
    memory so ``make_workdir`` can write the whole run, from its first line, to the step's log.

    Args:
        debug: Also show debug messages and tracebacks (``-v/--verbose``).
        quiet: Show only warnings and errors (``-q/--quiet``); ``debug`` wins if both are set.
    """
    global _log_buffer, _step_log, _step_log_used, _last_step_log
    file_level = logging.DEBUG if debug else logging.INFO
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG if debug else logging.WARNING if quiet else logging.INFO)
    console.setFormatter(_StatusFormatter(color=console.stream.isatty(), traceback=debug))
    _log_buffer = _BufferHandler(file_level)
    _step_log, _step_log_used, _last_step_log = None, False, None
    package_logger = logging.getLogger("aaftf")
    package_logger.handlers[:] = [console, _log_buffer]
    package_logger.setLevel(file_level)
    package_logger.propagate = False


def finish_logging(command: str | None, debug: bool) -> None:
    """Close the step log, if one is still open, and write out the messages still held.

    Messages logged after ``cleanup_workdir`` (e.g. the "next command" hint) are appended to the
    last step log if it still exists (the work directory was kept). A run whose subcommand never
    opened a step log (it has no work directory, e.g. ``trim``, ``sort``, ``assess``) writes the
    held messages to ``./<command>.log`` when ``debug`` is set.

    Args:
        command: Subcommand name, used for the fallback log file name.
        debug: Whether ``-v/--verbose`` was given.
    """
    _close_step_log(restore_buffer=False)
    target = None
    if _log_buffer is not None and _log_buffer.records:
        if _step_log_used:
            target = _last_step_log if _last_step_log is not None and _last_step_log.exists() else None
        elif debug and command:
            target = Path(f"{command}.log")
    if target is not None and _log_buffer is not None:
        handler = _plain_file_handler(target, _log_buffer.level)
        for record in _log_buffer.records:
            handler.handle(record)
        handler.close()
    if _log_buffer is not None:
        _log_buffer.records.clear()


def make_workdir(workdir: str | None, prefix: str) -> tuple[str, bool]:
    """Create a subcommand's working directory and start its log file there.

    When logging was set up by ``setup_logging`` (i.e. running from the CLI), the messages so far
    and everything logged until ``cleanup_workdir`` go to ``<workdir>/<prefix>.log`` (appended to,
    so steps sharing a ``--workdir`` each keep their own file across reruns).

    Args:
        workdir: User-supplied ``--workdir``, or None to auto-name it
            ``aaftf-<prefix>_<random id>``.
        prefix: Subcommand name used in the auto-generated name.

    Returns:
        Tuple (workdir, custom_workdir); pass both to ``cleanup_workdir``.
    """
    custom_workdir = bool(workdir)
    folder = workdir or f"aaftf-{prefix}_{uuid.uuid4().hex[:8]}"
    Path(folder).mkdir(parents=True, exist_ok=True)
    _open_step_log(Path(folder, f"{prefix}.log"))
    return folder, custom_workdir


def cleanup_workdir(workdir: str | Path, debug: bool, custom_workdir: bool) -> None:
    """Close the step log, then remove an auto-generated workdir unless debugging or the user supplied it.

    The step log lives in the workdir, so it is kept exactly when the workdir is.

    A user-supplied ``--workdir`` (which may be shared, e.g. by the ``pipeline_*`` subcommands, or
    even the current directory) is never deleted.

    Args:
        workdir: Working directory from ``make_workdir``.
        debug: Keep the directory for inspection.
        custom_workdir: Whether the user supplied ``workdir``.
    """
    _close_step_log()
    if not debug and not custom_workdir:
        safe_remove(workdir)


def close_step_log() -> None:
    """Close the open step log, if any, keeping the work directory and the log file.

    Call it from a ``finally`` after ``make_workdir`` so a step that raises still closes its log;
    after ``cleanup_workdir`` has run it does nothing.
    """
    _close_step_log()


class _Redirect308Handler(urllib.request.HTTPRedirectHandler):
    """Extend urllib's redirect handler to also follow HTTP 308.

    Python < 3.11 does not handle 308 (Permanent Redirect): the base
    ``redirect_request()`` only allows {301,302,303,307} and raises HTTPError
    for anything else, so both it and an http_error_308 dispatcher are needed.
    """

    def redirect_request(self, req: urllib.request.Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage, newurl: str) -> urllib.request.Request | None:
        """Treat a 308 as a 307 and build the redirected request with the base handler.

        Args:
            req: Original request.
            fp: Response body.
            code: HTTP status code.
            msg: HTTP reason phrase.
            headers: Response headers.
            newurl: Redirect target URL.

        Returns:
            The new request, or None if the base handler declines to redirect.
        """
        if code == 308:
            code = 307  # method-preserving permanent redirect
        return super().redirect_request(req, fp, code, msg, headers, newurl)

    def http_error_308(self, req: urllib.request.Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> Any:
        """Handle an HTTP 308 response the same way as a 302.

        Args:
            req: Original request.
            fp: Response body.
            code: HTTP status code.
            msg: HTTP reason phrase.
            headers: Response headers.

        Returns:
            The response for the redirected request.
        """
        return self.http_error_302(req, fp, code, msg, headers)


class _StatusFormatter(logging.Formatter):
    """Format records as ``[Mon DD HH:MM AM] message``; warnings and errors get a level prefix."""

    def __init__(self, color: bool, traceback: bool = True) -> None:
        """Set the date format, whether to colour the timestamp, and whether to show tracebacks.

        Args:
            color: Colour the timestamp green (for terminals).
            traceback: Append the traceback of records logged with ``exc_info``.
        """
        super().__init__(datefmt="%b %d %I:%M %p")
        self.color = color
        self.traceback = traceback

    def format(self, record: logging.LogRecord) -> str:
        """Format one log record.

        Args:
            record: The log record.

        Returns:
            The formatted message line.
        """
        stamp = f"[{self.formatTime(record, self.datefmt)}]"
        if self.color:
            stamp = f"\033[92m{stamp}\033[00m"
        message = record.getMessage()
        if record.levelno >= logging.WARNING:
            message = f"{record.levelname}: {message}"
        if self.traceback and record.exc_info:
            message = f"{message}\n{self.formatException(record.exc_info)}"
        return f"{stamp} {message}"


class _BufferHandler(logging.Handler):
    """Keep log records in memory until a step log file is opened to receive them."""

    def __init__(self, level: int) -> None:
        """Start with an empty buffer.

        Args:
            level: Lowest level to keep.
        """
        super().__init__(level)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        """Keep ``record``.

        Args:
            record: The log record.
        """
        self.records.append(record)


def _plain_file_handler(path: Path, level: int) -> logging.FileHandler:
    """Return an appending log file handler without colour, including tracebacks.

    Args:
        path: Log file path.
        level: Lowest level written.

    Returns:
        The file handler.
    """
    handler = logging.FileHandler(path, mode="a")
    handler.setLevel(level)
    handler.setFormatter(_StatusFormatter(color=False, traceback=True))
    return handler


def _open_step_log(path: Path) -> None:
    """Write the held messages to ``path`` and send further messages there too.

    Does nothing unless ``setup_logging`` ran (e.g. when a ``run()`` is called directly from Python).

    Args:
        path: The step log file, inside the work directory.
    """
    global _step_log, _step_log_used, _last_step_log
    if _log_buffer is None:
        return
    _close_step_log(restore_buffer=False)
    handler = _plain_file_handler(path, _log_buffer.level)
    for record in _log_buffer.records:
        handler.handle(record)
    _log_buffer.records.clear()
    package_logger = logging.getLogger("aaftf")
    package_logger.removeHandler(_log_buffer)
    package_logger.addHandler(handler)
    _step_log, _step_log_used, _last_step_log = handler, True, path


def _close_step_log(restore_buffer: bool = True) -> None:
    """Close the open step log, if any, and go back to holding messages in memory.

    Args:
        restore_buffer: Re-attach the memory buffer so messages before the next step's log (e.g.
            the next pipeline step) are not lost.
    """
    global _step_log
    package_logger = logging.getLogger("aaftf")
    if _step_log is not None:
        package_logger.removeHandler(_step_log)
        _step_log.close()
        _step_log = None
    if restore_buffer and _log_buffer is not None and _log_buffer not in package_logger.handlers:
        package_logger.addHandler(_log_buffer)


def _count_lines(fh: BinaryIO | IO[bytes]) -> int:
    """Count lines in a binary stream, including a final line with no newline.

    Args:
        fh: Binary stream read in 1 MiB chunks until EOF.

    Returns:
        Number of lines.
    """
    lines = 0
    last = b"\n"
    for chunk in iter(lambda: fh.read(1 << 20), b""):
        lines += chunk.count(b"\n")
        last = chunk[-1:]
    return lines + (last != b"\n")


@functools.cache
def _url_opener() -> urllib.request.OpenerDirector:
    """Build (once) the urllib opener used by ``open_url``.

    Returns:
        An opener that follows HTTP 308 redirects.
    """
    return urllib.request.build_opener(_Redirect308Handler())
