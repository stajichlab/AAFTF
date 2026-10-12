"""Remove contaminant contigs from one or more assemblies with NCBI FCS-GX.

FCS-GX reads its whole database (~470 GB for the full gxdb) in a random-access pattern, so it
needs that much RAM and fast storage; this is a standalone utility for large-memory nodes rather
than a pipeline step. A database on tmpfs (``/dev/shm``) or on storage that passes a quick
random-read test is used in place; otherwise it is copied once into a per-job folder under
``$TMPDIR`` (``--copy_to``), every assembly is screened against it, and the copy is removed at
the end.
"""

import logging
import os
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

from aaftf.utility import available_cpus, cleanup_workdir, close_step_log, fasta_stats, filter_fasta, get_ram, make_workdir, print_cmd

__all__ = ["run", "cleaned_name", "read_latency_ms"]


logger = logging.getLogger(__name__)

# parent folder of the per-job database copy when --copy_to is not given and $TMPDIR is not set
FALLBACK_COPY_TO = "/scratch"

# memory needed on top of the database for run_gx.py itself, as a fraction of the database size
MEMORY_HEADROOM = 0.05

# random-read speed test: a database whose mean 4 KiB random read takes longer than this is copied
FAST_READ_MS = 0.2
SPEEDTEST_READS = 256
SPEEDTEST_MAX_SECONDS = 3.0

# filesystems reached over the network: a database copy there is no faster than the original
NETWORK_FS = ("lustre", "gpfs", "nfs", "nfs4", "beegfs", "cifs", "smb3", "fuse.sshfs", "ceph", "panfs", "wekafs")

_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas", ".fsa")


def _prefix(path: str) -> str:
    """Return the file name of ``path`` without its FASTA (and ``.gz``) extension."""
    name = Path(path).name
    name = name.removesuffix(".gz")
    for suffix in _FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def cleaned_name(path: str, outdir: str | None = None) -> str:
    """Return the cleaned-assembly path for ``path``: ``<prefix>.fcsgx_cleaned.fasta`` in ``outdir`` (default: next to ``path``)."""
    return str(Path(outdir or Path(path).parent, f"{_prefix(path)}.fcsgx_cleaned.fasta"))


def _resolve_db(db: str) -> str:
    """Return the FCS-GX database path prefix for ``db``: a prefix (``gxdb/all``), its ``.gxi`` file, or the folder holding it.

    Raises:
        FileNotFoundError: If no ``.gxi`` is found, or a folder holds several.
    """
    path = Path(db)
    if path.is_dir():
        gxi = sorted(path.glob("*.gxi"))
        if len(gxi) != 1:
            found = ", ".join(f.name for f in gxi) or "none"
            raise FileNotFoundError(f"-d/--db folder {db} must hold exactly one FCS-GX .gxi index (found: {found}); give the database path prefix instead (e.g. {path}/all)")
        return str(gxi[0].with_suffix(""))
    db = db.removesuffix(".gxi")
    if not Path(f"{db}.gxi").is_file():
        raise FileNotFoundError(f"FCS-GX database {db}.gxi not found; give the database path prefix with -d/--db (e.g. /path/gxdb/all)")
    return db


def _db_files(db: str) -> list[Path]:
    """Return every file of the FCS-GX database with path prefix ``db`` (``{db}.gxi``, ``{db}.gxs``, ``{db}.taxa.tsv``, ...)."""
    prefix = Path(db)
    return sorted(f for f in prefix.parent.glob(f"{prefix.name}.*") if f.is_file())


def _fs_type(path: str) -> str:
    """Return the filesystem type (``tmpfs``, ``xfs``, ``lustre``, ...) holding ``path``, or "" if unknown."""
    target = os.path.realpath(path)
    best, fstype = "", ""
    try:
        with open("/proc/mounts") as fh:
            for line in fh:
                fields = line.split()
                if len(fields) < 3:
                    continue
                mount = fields[1].replace("\\040", " ")
                if (target == mount or target.startswith(mount.rstrip("/") + "/")) and len(mount) > len(best):
                    best, fstype = mount, fields[2]
    except OSError:
        pass
    return fstype


def _in_ram(path: str) -> bool:
    """Return whether ``path`` is on a RAM-backed filesystem (tmpfs, e.g. ``/dev/shm``)."""
    return _fs_type(path) in ("tmpfs", "ramfs")


def read_latency_ms(path: str, reads: int = SPEEDTEST_READS, max_seconds: float = SPEEDTEST_MAX_SECONDS, block: int = 4096) -> float:
    """Return the mean time in ms of a random ``block``-byte read from ``path``, the access pattern of FCS-GX.

    Stops after ``reads`` reads or ``max_seconds``, whichever comes first, so slow storage is not
    tested for long. Data already in the page cache is (correctly) measured as fast.
    """
    size = os.path.getsize(path)
    blocks = max(size // block, 1)
    rng = random.Random()
    fd = os.open(path, os.O_RDONLY)
    try:
        done, start = 0, time.perf_counter()
        while done < reads:
            os.pread(fd, block, rng.randrange(blocks) * block)
            done += 1
            if time.perf_counter() - start > max_seconds:
                break
        elapsed = time.perf_counter() - start
    finally:
        os.close(fd)
    return elapsed / done * 1000


def _needs_copy(db: str, files: list[Path]) -> bool:
    """Return whether to copy the database: not when it is on tmpfs or its largest file passes the random-read speed test."""
    if _in_ram(db):
        logger.info(f"FCS-GX database {db} is on tmpfs; using it in place")
        return False
    largest = max(files, key=lambda f: f.stat().st_size)
    latency = read_latency_ms(str(largest))
    fstype = _fs_type(db) or "unknown filesystem"
    if latency <= FAST_READ_MS:
        logger.info(f"FCS-GX database {db} ({fstype}) reads in {latency:.3f} ms per random 4 KiB block (<= {FAST_READ_MS} ms); using it in place")
        return False
    logger.info(f"FCS-GX database {db} ({fstype}) reads in {latency:.3f} ms per random 4 KiB block (> {FAST_READ_MS} ms); copying it")
    return True


def _same_filesystem(path: str, folder: str) -> bool:
    """Return whether ``folder`` (or, if it does not exist yet, its nearest existing parent) is on the same filesystem as ``path``."""
    existing = Path(folder).absolute()
    while not existing.exists() and existing != existing.parent:
        existing = existing.parent
    return os.stat(path).st_dev == os.stat(existing).st_dev


def _copy_target_ok(db_file: str, copy_to: str) -> bool:
    """Return whether copying the database into ``copy_to`` can help; warn when it cannot, or may not.

    A copy on the database's own filesystem cannot be faster, so it is skipped; a copy on a
    network filesystem is made but warned about.
    """
    if _same_filesystem(db_file, copy_to):
        logger.warning(f"--copy_to {copy_to} is on the same filesystem as the FCS-GX database, so a copy would not be faster; using the database in place. Give a node-local --copy_to to copy it")
        return False
    fstype = _fs_type(copy_to)
    if fstype in NETWORK_FS:
        logger.warning(f"--copy_to {copy_to} is on a network filesystem ({fstype}); the database copy may be no faster than the original. Prefer node-local storage")
    return True


def _check_memory(db_bytes: int, in_ram: bool) -> None:
    """Raise if this job does not have the RAM to keep a ``db_bytes`` database cached plus run_gx.py's own memory.

    ``get_ram()`` covers SLURM / cgroup limits. A database already on tmpfs is resident, so it is
    not checked.

    Raises:
        MemoryError: If the RAM available to this job is too small.
    """
    gib = 1024.0**3
    if in_ram:
        return
    need_gb = db_bytes * (1 + MEMORY_HEADROOM) / gib
    avail_gb = get_ram()
    logger.info(f"FCS-GX database is {db_bytes / gib:,.1f} GB; needs ~{need_gb:,.0f} GB RAM, {avail_gb:,.0f} GB available")
    if avail_gb < need_gb:
        raise MemoryError(f"FCS-GX needs ~{need_gb:,.0f} GB of RAM for its database, but only {avail_gb:,.0f} GB is available to this job. Run on a large-memory node (e.g. sbatch --mem={int(need_gb) + 10}G).")


def _copy_db(db: str, files: list[Path], copy_to: str) -> tuple[str, str]:
    """Copy the database ``files`` into a new per-job folder in ``copy_to``; return (that folder, the database prefix in it).

    Raises:
        OSError: If ``copy_to`` does not have room for the database.
    """
    db_bytes = sum(f.stat().st_size for f in files)
    Path(copy_to).mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(copy_to).free
    if free < db_bytes:
        raise OSError(f"{copy_to} has only {free / 1024**3:,.0f} GB free, less than the {db_bytes / 1024**3:,.0f} GB FCS-GX database; choose another --copy_to")
    job = os.environ.get("SLURM_JOB_ID", str(os.getpid()))
    dest = tempfile.mkdtemp(prefix=f"aaftf-gxdb-{job}-", dir=copy_to)
    logger.info(f"Copying the FCS-GX database to {dest}")
    start = time.monotonic()
    for f in files:
        shutil.copy2(f, dest)
    logger.info(f"Copied the FCS-GX database in {time.monotonic() - start:,.0f} s")
    return dest, str(Path(dest, Path(db).name))


def _screen(assembly: str, outfile: str, report: str, run_dir: Path, taxid: int, db: str, cpus: int, debug: bool) -> None:
    """Run run_gx.py on ``assembly`` and write the contigs it does not mark ``EXCLUDE`` to ``outfile``.

    Other actions (``TRIM``, ``FIX``, ``REVIEW``, ...) are only logged as warnings, since removing
    part of a contig is not supported here. The FCS-GX report is copied to ``report``.

    Raises:
        RuntimeError: If ``run_gx.py`` fails or writes no report.
    """
    run_dir.mkdir(parents=True, exist_ok=True)
    num_seqs, assembly_size = fasta_stats(assembly)
    logger.info(f"{assembly}: {num_seqs:,} contigs and {assembly_size:,} bp")

    fcsgx_compute = ["run_gx.py", "--fasta", assembly, "--tax-id", f"{taxid}", "--gx-db", db, "--out-dir", str(run_dir)]
    print_cmd(fcsgx_compute)
    fcs_log = run_dir / "fcs_gx.log"
    with open(fcs_log, "w") as logfile:
        result = subprocess.run(fcsgx_compute, stderr=logfile, env={**os.environ, "GX_NUM_CORES": str(cpus)})

    # output tsv:
    # seq_id	start_pos	end_pos	seq_len	action	div	agg_cont_cov	top_tax_name
    reports = sorted(run_dir.glob("*.fcs_gx_report.txt"))
    if result.returncode != 0 or not reports:
        raise RuntimeError(f"run_gx.py failed (exit {result.returncode}) or did not write a *.fcs_gx_report.txt for {assembly}; see {fcs_log}")
    fcsgx_tsv = reports[0]

    seqs_to_drop = set()
    with open(fcsgx_tsv) as fcsgx_out:
        for line in fcsgx_out:
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(cols) < 5:
                continue
            action = cols[4].removeprefix("ACTION_")
            if action == "EXCLUDE":
                seqs_to_drop.add(cols[0])
            else:
                logger.warning(f"FCS-GX reports {action} for {cols[0]}:{cols[1]}..{cols[2]}; not applied, review it manually")

    logger.info(f"Dropping {len(seqs_to_drop)} contigs flagged by FCS-GX")
    num_seqs, assembly_size = filter_fasta(assembly, outfile, lambda seq_id: seq_id not in seqs_to_drop)
    if debug:
        print("Contigs dropped by FCS-GX: {:}".format(",".join(sorted(seqs_to_drop))))
    logger.info(f"{outfile}: {num_seqs:,} contigs and {assembly_size:,} bp")
    shutil.copy(fcsgx_tsv, report)


def run(
    input: list[str],
    db: str,
    taxid: int,
    outdir: str | None = None,
    workdir: str | None = None,
    cpus: int | None = None,
    copy_to: str | None = None,
    debug: bool = False,
    **kwargs: Any,
) -> None:
    """Run NCBI FCS-GX on each assembly in ``input`` and drop the contigs it marks ``EXCLUDE``.

    Before anything else, checks that the job has enough RAM (``get_ram()``, which respects SLURM
    limits) to keep the database cached. The database is used in place when it is on tmpfs (e.g.
    ``/dev/shm``) or passes a random-read speed test; otherwise it is copied once into a per-job
    folder in ``copy_to`` and removed at the end, also when a run fails. For
    each input, writes ``<prefix>.fcsgx_cleaned.fasta`` and the FCS-GX report
    ``<prefix>.fcs_gx_report.txt``.

    Args:
        input: Assembly FASTA files to screen.
        db: FCS-GX database: its path prefix (``gxdb/all``), its ``.gxi`` file, or the folder holding it.
            Every file sharing the prefix is part of the database.
        taxid: NCBI taxonomy ID of the expected organism (or a parent taxon).
        outdir: Directory for the outputs; None writes each next to its input.
        workdir: Working directory; a temporary one is created when None.
        cpus: Threads for run_gx.py (``GX_NUM_CORES``); None uses ``available_cpus()``.
        copy_to: Folder for the per-job database copy (fast node-local storage); None uses
            ``$TMPDIR``, or ``FALLBACK_COPY_TO`` when that is not set.
        debug: Print dropped contigs and keep the work directory.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.

    Raises:
        FileNotFoundError: If ``run_gx.py`` is not on PATH, or an input or the database's ``.gxi`` does not exist.
        ValueError: If two inputs would write the same output file.
        MemoryError: If there is not enough RAM for the database.
        OSError: If ``copy_to`` has no room for the database copy.
        RuntimeError: If ``run_gx.py`` fails or writes no report.
    """
    if shutil.which("run_gx.py") is None:
        raise FileNotFoundError("run_gx.py (NCBI FCS-GX) not found on PATH; install FCS-GX or load its environment (see AAFTF dependency)")
    if isinstance(input, str):
        input = [input]
    missing = [f for f in input if not Path(f).is_file()]
    if missing:
        raise FileNotFoundError(f"Input assembly not found: {', '.join(missing)}")
    outputs = [cleaned_name(f, outdir) for f in input]
    if len(set(outputs)) != len(outputs):
        raise ValueError("Two inputs would write the same .fcsgx_cleaned.fasta; give them different names or run them separately")
    db = _resolve_db(db)
    db_files = _db_files(db)
    copy_to = copy_to or os.environ.get("TMPDIR") or FALLBACK_COPY_TO
    cpus = cpus or available_cpus()

    workdir, custom_workdir = make_workdir(workdir, "fcs_gx_purge")
    db_copy = None
    try:
        _check_memory(sum(f.stat().st_size for f in db_files), _in_ram(db))
        gx_db = db
        if _needs_copy(db, db_files) and _copy_target_ok(str(db_files[0]), copy_to):
            db_copy, gx_db = _copy_db(db, db_files, copy_to)

        if outdir:
            Path(outdir).mkdir(parents=True, exist_ok=True)
        for n, (assembly, outfile) in enumerate(zip(input, outputs, strict=True), 1):
            logger.info(f"Screening assembly {n} of {len(input)}: {assembly}")
            report = str(Path(outfile).parent / f"{_prefix(assembly)}.fcs_gx_report.txt")
            _screen(assembly, outfile, report, Path(workdir, f"{n}_{_prefix(assembly)}"), taxid, gx_db, cpus, debug)

        cleanup_workdir(workdir, debug, custom_workdir)
    finally:
        if db_copy:
            logger.info(f"Removing the FCS-GX database copy {db_copy}")
            shutil.rmtree(db_copy, ignore_errors=True)
        close_step_log()
