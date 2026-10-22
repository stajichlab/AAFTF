"""Download and cache AAFTF reference databases.

``AAFTF database`` with no arguments lists every database in
``aaftf.resources.DATABASES``: its abbreviation, file, size, which subcommands
use it, and where it is stored. Naming databases (by abbreviation or file name,
or ``all``) downloads them into the database folder (see
``utility.db_dirs``), so later AAFTF commands reuse the cached files.
"""

import contextlib
import logging
import os
import tarfile
import urllib.request
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from Bio import SeqIO

from aaftf.resources import CONTAMDB, CONTAMDB_BBTOOLS, DATABASES, TAXPURGE_DEFAULT_DATABASES, VECSCREEN_EXCLUDE_TERMS
from aaftf.utility import bbtools_resources, db_dirs, db_file, db_write_dir, download_file, find_db_file, find_db_index, open_maybe_gz, open_url, require_tools, run_cmd, safe_remove, warn_if_home_cache, write_fasta

__all__ = ["build_contamdb", "find_contamdb", "require_contamdb", "run"]

# files makeblastdb writes for a (single-volume) nucleotide database
_BLASTDB_SUFFIXES = (".nhr", ".nin", ".nsq")


logger = logging.getLogger(__name__)


def run(databases: list[str] | None = None, force: bool = False, no_build: bool = False, **kwargs: Any) -> None:
    """Execute the ``database`` subcommand.

    Each sourmash database is prepared as a RocksDB index, its final form: a published index is
    downloaded and unpacked, otherwise the sketch file is downloaded, indexed, and deleted.

    Args:
        databases: Database abbreviations or file names (case-insensitive), or ``required``
            (what the default pipeline needs), ``optional`` (the rest), ``taxpurge`` (its default
            sourmash databases) or ``all``. When empty, list the databases instead of downloading.
        force: Re-download (and rebuild) even if a copy already exists (into the first
            writable database folder).
        no_build: For a sourmash database without a published index, download its sketch file
            only, without building (and replacing it with) the RocksDB index.

    A file that was only needed to build a database that is already built (a sourmash sketch
    file once its index exists, or PhiX/UniVec_Core once ``contamdb.fa`` exists) is skipped
    unless ``force`` is given.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``debug``, ...); ignored.

    Raises:
        ValueError: If a database name is not recognised.
        RuntimeError: If any download or index build fails.
    """
    if not databases:
        _list_db()
        return

    selected = _resolve(databases)
    db_dir = db_write_dir()
    warn_if_home_cache()
    logger.info(f"AAFTF database directory: {db_dir}")

    errors = []
    for abbr in selected:
        entry = DATABASES[abbr]
        try:
            fulfilled = None if force else _fulfilled_by(abbr)
            if fulfilled:
                logger.info(f"{abbr}: fulfilled by {fulfilled}, skipped (--force downloads it)")
            elif entry.get("built_from"):
                build_contamdb(force)
            elif entry.get("index"):
                _install_index_archive(entry, force)
            elif entry.get("build_index") and not no_build:
                _build_index(entry, force)
            else:
                dest = db_file(entry["filename"], force=force)
                already_present = Path(dest).exists() and not force
                download_file(entry["url"], dest, force=force)
                if entry.get("executable") and not already_present:
                    os.chmod(dest, 0o555)
        except Exception as e:
            errors.append(f"{abbr} ({entry['filename']}): {e}")

    if errors:
        hint = "Re-run the same command later; files already downloaded are skipped unless you pass --force."
        raise RuntimeError("Some downloads failed:\n" + "\n".join(errors) + f"\n{hint}")
    logger.info("Download complete. Run 'AAFTF database' to see where each database is stored.")


def build_contamdb(force: bool = False) -> str:
    """Build ``contamdb.fa`` and its BLAST database, the contamination database of ``filter`` and ``vecscreen``.

    The FASTA is UniVec_Core without the entries whose names contain any of
    ``VECSCREEN_EXCLUDE_TERMS`` (fungal genomic sequence that would match, and split, real fungal
    contigs), followed by the PhiX genome and the BBTools ``lambda.fa.gz`` and ``adapters.fa``.
    UniVec_Core and PhiX are downloaded first if needed. The FASTA is written to a temporary file
    and renamed when complete, then indexed by ``require_contamdb(blast=True)``. An existing
    ``contamdb.fa`` (which users may have added their own sequences to) is kept; only its BLAST
    database is brought up to date. Once the database is built, UniVec_Core and PhiX are deleted,
    but only from the folder AAFTF downloads into (a copy in another ``$AAFTF_DB`` folder may be
    shared).

    Args:
        force: Rebuild the FASTA (and re-download UniVec_Core and PhiX) even if it is present.

    Returns:
        Path prefix of the database (``<prefix>.fa`` is the FASTA; ``blastn -db <prefix>``).

    Raises:
        FileNotFoundError: If ``makeblastdb``, ``bbduk.sh`` or a BBTools resource file is missing.
        RuntimeError: If ``makeblastdb`` fails or does not produce a complete database.
    """
    existing = find_contamdb()
    if existing and not force:
        logger.info(f"Already present: {existing}.fa")
        return require_contamdb(blast=True)
    univec, phix = (download_file(DATABASES[abbr]["url"], db_file(DATABASES[abbr]["filename"], force=force), force=force) for abbr in DATABASES["contamdb"]["built_from"])
    bbtools = bbtools_resources(CONTAMDB_BBTOOLS)

    fasta = db_write_dir() / CONTAMDB["fasta"]
    tmp = fasta.with_name(fasta.name + ".tmp")
    kept = dropped = 0
    with open(tmp, "w") as out:
        with open_maybe_gz(univec) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                if any(term in record.description.lower() for term in VECSCREEN_EXCLUDE_TERMS):
                    dropped += 1
                    continue
                write_fasta(out, record.description, str(record.seq))
                kept += 1
        for source in (phix, *bbtools):
            with open_maybe_gz(source) as fh:
                for record in SeqIO.parse(fh, "fasta"):
                    write_fasta(out, record.description, str(record.seq))
    tmp.replace(fasta)
    logger.info(f"Kept {kept:,} UniVec_Core sequences ({dropped:,} fungal/yeast-derived left out), added PhiX, {', '.join(CONTAMDB_BBTOOLS)}: {fasta}")
    prefix = require_contamdb(blast=True)
    for source in (univec, phix):
        if Path(source).parent == db_write_dir():
            safe_remove(source)
            logger.info(f"Removed {source}; {CONTAMDB['fasta']} replaces it")
    return prefix


def find_contamdb() -> str | None:
    """Return the path prefix of the first ``contamdb.fa`` in ``db_dirs()`` (without ``.fa``), or None."""
    for folder in db_dirs():
        if (folder / CONTAMDB["fasta"]).is_file():
            return str(folder / CONTAMDB["blastdb"])
    return None


def require_contamdb(blast: bool = False) -> str:
    """Return the path prefix of the contamination database, raising if ``contamdb.fa`` is missing.

    Only ``contamdb.fa`` is required, so users can add their own sequences to it. With ``blast``,
    its BLAST database (``.nhr``, ``.nin``, ``.nsq`` next to it) is also checked, and rebuilt with
    ``makeblastdb`` when a file is missing, empty or older than ``contamdb.fa``.

    Args:
        blast: Also make sure the BLAST database is complete and current.

    Returns:
        The path prefix: ``<prefix>.fa`` is the FASTA, and ``blastn -db <prefix>`` the BLAST database.

    Raises:
        FileNotFoundError: If no ``contamdb.fa`` is in any ``db_dirs()`` folder, or ``makeblastdb`` is missing.
        RuntimeError: If ``makeblastdb`` fails or does not produce a complete database.
    """
    prefix = find_contamdb()
    if prefix is None:
        searched = ", ".join(str(folder) for folder in db_dirs())
        raise FileNotFoundError(f"missing database: contamdb ({CONTAMDB['fasta']}; searched {searched})\nBuild it first: AAFTF database contamdb")
    if blast and not _blastdb_current(prefix):
        require_tools(["makeblastdb"], hint="NCBI BLAST+ builds the contamdb BLAST database")
        result = run_cmd(["makeblastdb", "-dbtype", "nucl", "-in", prefix + ".fa", "-out", prefix, "-title", CONTAMDB["blastdb"]])
        if result.returncode != 0 or not _blastdb_current(prefix):
            raise RuntimeError(f"makeblastdb failed (exit status {result.returncode}) building the BLAST database {prefix}")
    return prefix


def _blastdb_current(prefix: str) -> bool:
    """Return True if every BLAST database file of ``prefix`` is present, non-empty and not older than ``<prefix>.fa``."""
    fasta_mtime = Path(prefix + ".fa").stat().st_mtime
    files = [Path(prefix + suffix) for suffix in _BLASTDB_SUFFIXES]
    return all(f.is_file() and f.stat().st_size > 0 and f.stat().st_mtime >= fasta_mtime for f in files)


def _fulfilled_by(abbr: str) -> str | None:
    """Return what makes database ``abbr``'s own file unnecessary, or None.

    That is ``"rocksdb"`` for a sourmash sketch file whose RocksDB index is built, or the name of
    a locally built database (e.g. ``"contamdb"``) that lists ``abbr`` in ``built_from`` and is
    built. A file still present in a database folder is not reported, since it is in use.

    Args:
        abbr: Key of ``DATABASES``.

    Returns:
        ``"rocksdb"``, the building database's abbreviation, or None.
    """
    entry = DATABASES[abbr]
    if entry.get("index") or find_db_file(entry["filename"]):
        return None
    if entry.get("build_index") and find_db_index(entry["build_index"]):
        return "rocksdb"
    return next((name for name, other in DATABASES.items() if abbr in other.get("built_from", []) and find_db_file(other["filename"])), None)


def _install_index_archive(entry: dict[str, Any], force: bool) -> None:
    """Download a published RocksDB index archive and unpack it into the database folder.

    The archive is downloaded and unpacked in a temporary folder, and the index moved into place
    only once complete (it must hold RocksDB's ``CURRENT`` file), so an interrupted run never
    leaves a partial index that a later run, or ``taxpurge``, would use. The archive is deleted.

    Args:
        entry: The ``DATABASES`` entry: ``filename`` is the index folder, ``url`` its archive.
        force: Re-install even if the index is already present.

    Raises:
        RuntimeError: If the archive does not hold the expected index folder.
    """
    existing = find_db_index(entry["filename"])
    if existing and not force:
        logger.info(f"Already present: {existing}")
        return
    with _staging(entry["filename"]) as tmp:
        archive = download_file(entry["url"], str(tmp / Path(entry["url"]).name))
        logger.info(f"Unpacking {Path(archive).name} ...")
        _extract_archive(archive, tmp)
        safe_remove(archive)
        _move_index(tmp / entry["filename"])


def _build_index(entry: dict[str, Any], force: bool) -> None:
    """Build a sourmash database's RocksDB index from its sketch file, then delete the sketch file.

    The sketch file is downloaded first if needed. The index is built in a temporary folder and
    moved into place only once complete. The sketch file is deleted only from the folder AAFTF
    downloads into (a copy in another ``$AAFTF_DB`` folder may be shared).

    Args:
        entry: The ``DATABASES`` entry: ``filename`` / ``url`` are the sketch file, ``build_index``
            the index folder, ``scaled`` its resolution.
        force: Rebuild (and re-download the sketch file) even if the index is already present.

    Raises:
        RuntimeError: If the build does not produce a complete index.
    """
    existing = find_db_index(entry["build_index"])
    if existing and not force:
        logger.info(f"Already present: {existing}")
        return
    sketches = download_file(entry["url"], db_file(entry["filename"], force=force), force=force)
    require_tools(["sourmash"], hint="sourmash and its branchwater plugin (sourmash_plugin_branchwater) build the index; or use --no-build to keep the sketch file")
    with _staging(entry["build_index"]) as tmp:
        logger.info(f"Building {entry['build_index']} from {Path(sketches).name} (this can take a while)")
        run_cmd(["sourmash", "scripts", "index", sketches, "-o", str(tmp / entry["build_index"]), "-k", "31", "-s", str(entry["scaled"])], False)
        _move_index(tmp / entry["build_index"])
    if Path(sketches).parent == db_write_dir():
        safe_remove(sketches)
        logger.info(f"Removed {sketches}; the RocksDB index replaces it")


@contextlib.contextmanager
def _staging(name: str) -> Iterator[Path]:
    """Yield an empty temporary folder next to where index ``name`` goes, removed afterwards."""
    tmp = db_write_dir() / f"{name}.tmp"
    safe_remove(tmp)
    tmp.mkdir(parents=True)
    try:
        yield tmp
    finally:
        safe_remove(tmp)


def _move_index(built: Path) -> None:
    """Move a complete RocksDB index ``built`` into the database folder (replacing an old copy).

    Raises:
        RuntimeError: If ``built`` is not a complete index (no ``CURRENT`` file).
    """
    if not (built / "CURRENT").is_file():
        raise RuntimeError(f"no complete RocksDB index {built.name} was produced")
    dest = db_write_dir() / built.name
    safe_remove(dest)
    built.rename(dest)
    logger.info(f"Saved {dest}")


def _extract_archive(archive: str, dest: Path) -> None:
    """Extract tar archive ``archive`` into ``dest``, refusing members that would land outside it.

    Uses tarfile's ``data`` filter where available (Python 3.12, and 3.10/3.11 patch releases);
    otherwise rejects absolute paths, ``..`` components and links before extracting.

    Raises:
        RuntimeError: If the archive has an unsafe member.
    """
    with tarfile.open(archive) as tar:
        if hasattr(tarfile, "data_filter"):
            tar.extractall(dest, filter="data")
            return
        for member in tar.getmembers():
            name = Path(member.name)
            if name.is_absolute() or ".." in name.parts or member.issym() or member.islnk():
                raise RuntimeError(f"unsafe path in {Path(archive).name}: {member.name}")
        tar.extractall(dest)


def _folder_size(path: str) -> int:
    """Return the total size of the files in folder ``path`` (a RocksDB index; not recursive)."""
    return sum(f.stat().st_size for f in Path(path).iterdir() if f.is_file())


def _list_db() -> None:
    """List every database with its abbreviation, file, size, users and storage folder.

    Databases are listed in two groups, like ``AAFTF dependency``: those ``AAFTF pipeline_short`` needs
    with its default settings, and those only needed by non-default options or optional steps.
    Missing databases show their remote (estimated) size and "not downloaded".
    """
    folders = db_dirs()
    logger.info("Database folders, searched in order:" if len(folders) > 1 else "Database folder:")
    for folder in folders:
        logger.info(f"{folder}{'' if folder.is_dir() else ' (does not exist yet)'}")

    rows: list[tuple[str, str, str, str, str]] = []  # (name, file, size label, used by, location)
    section_starts = {}  # row index -> heading printed before it
    downloaded_total = 0
    pending_total = 0
    pending_unknown = False
    for required, heading in ((True, "Needed by the default pipeline:"), (False, "Only for non-default options or optional steps:")):
        section_starts[len(rows)] = heading
        for abbr, entry in DATABASES.items():
            if entry["required"] is not required:
                continue
            index = find_db_index(entry["build_index"]) if entry.get("build_index") else None
            local_path = find_db_index(entry["filename"]) if entry.get("index") else find_db_file(entry["filename"])
            size: int | None
            if local_path:
                size = _folder_size(local_path) if entry.get("index") else Path(local_path).stat().st_size
                downloaded_total += size
                size_label, location = _human_size(size), str(Path(local_path).parent)
            elif entry.get("built_from"):  # built locally by AAFTF database, nothing to download
                size_label, location = "-", "not built"
            elif fulfilled := _fulfilled_by(abbr):  # only needed to build an index/database that is built; its download size, not counted
                size = _remote_size(entry["url"])
                size_label, location = ("unknown" if size is None else _human_size(size)), f"fulfilled by {fulfilled}"
            else:
                size = _remote_size(entry["url"])
                if size is None:
                    pending_unknown = True
                else:
                    pending_total += size
                size_label, location = ("unknown" if size is None else _human_size(size)), "not downloaded"
            rows.append((abbr, entry["filename"], size_label, entry["used_by"], location))
            if entry.get("build_index"):  # the index AAFTF database builds from it, right below it
                if index:
                    size = _folder_size(index)
                    downloaded_total += size
                    rows.append(("", entry["build_index"], _human_size(size), "", str(Path(index).parent)))
                else:
                    rows.append(("", entry["build_index"], "-", "", "not built"))

    # Anything else directly in the database folders (e.g. leftover reports); subfolders are listed
    # as one entry, without their size and without listing their contents
    known_files = {entry["filename"] for entry in DATABASES.values()} | {entry["build_index"] for entry in DATABASES.values() if entry.get("build_index")}
    other_start = len(rows)
    for folder in folders:
        if not folder.is_dir():
            continue
        for f in sorted(f for f in folder.iterdir() if f.name not in known_files and not f.name.startswith(CONTAMDB["blastdb"] + ".n")):
            if f.is_dir():
                rows.append(("-", f"{f.name}/", "-", "-", f"{folder} (other folder)"))
            elif f.is_file():
                size = f.stat().st_size
                downloaded_total += size
                rows.append(("-", f.name, _human_size(size), "-", f"{folder} (other file)"))

    if len(rows) > other_start:
        section_starts[other_start] = "Other files in the database folders:"

    header = ("Name", "File", "Size", "Used by", "Location")
    name_w, file_w, size_w, used_w = (max(len(row[i]) for row in [header, *rows]) for i in range(4))
    print(f"  {header[0]:<{name_w}}  {header[1]:<{file_w}}  {header[2]:>{size_w}}  {header[3]:<{used_w}}  {header[4]}")
    for i, (name, filename, size_label, used_by, location) in enumerate(rows):
        if i in section_starts:
            print(section_starts[i])
        print(f"  {name:<{name_w}}  {filename:<{file_w}}  {size_label:>{size_w}}  {used_by:<{used_w}}  {location}")

    pending_label = _human_size(pending_total) + ("+" if pending_unknown else "")
    print("-" * (name_w + file_w + size_w + used_w + 20))
    print(f"  {'Downloaded':<32} {_human_size(downloaded_total):>10}")
    print(f"  {'Not yet downloaded (estimated)':<32} {pending_label:>10}")
    print("\nDownload with: AAFTF database NAME [NAME ...]   (abbreviation or file name, or 'required', 'optional', 'taxpurge', 'all')")
    print("Sourmash databases are prepared as RocksDB indexes (the sketch file is then removed); add --no-build to keep the sketch file instead.")


def _human_size(num_bytes: float) -> str:
    """Format a byte count as a human-readable string (e.g. ``1.2 GB``).

    Args:
        num_bytes: Size in bytes.

    Returns:
        The size in B/KB/MB/GB/TB (1024-based) with one decimal place.
    """
    size = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"  # pragma: no cover - unreachable, but keeps a return for every path


def _remote_size(url: str) -> int | None:
    """Return the remote file size in bytes via an HTTP HEAD request, or ``None`` if unavailable.

    Args:
        url: URL of the remote file.

    Returns:
        The ``Content-Length`` in bytes, or None if the header is missing or the request fails.
    """
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "AAFTF/1.0"})
        with open_url(req, timeout=10) as response:
            length = response.headers.get("Content-Length")
            return int(length) if length is not None else None
    except Exception:
        return None


def _resolve(names: list[str]) -> list[str]:
    """Map database abbreviations / file names (or a group keyword) to abbreviations, in order, without duplicates.

    Args:
        names: Database abbreviations or file names (case-insensitive), or ``required`` (what the
            default pipeline needs), ``optional`` (the rest), ``taxpurge`` (``taxpurge``'s default
            sourmash databases) or ``all``. A database's companions (``with``, e.g. a sourmash
            database's lineage CSV) are added after it.

    Returns:
        Keys of ``DATABASES`` in the order first requested.

    Raises:
        ValueError: If any name matches no database.
    """
    lookup = {}
    for abbr, entry in DATABASES.items():
        lookup[abbr.lower()] = abbr
        lookup[entry["filename"].lower()] = abbr

    selected: list[str] = []
    unknown: list[str] = []
    for name in names:
        key = name.lower()
        if key == "all":
            matches = list(DATABASES)
        elif key == "taxpurge":
            matches = list(TAXPURGE_DEFAULT_DATABASES)
        elif key in ("required", "optional"):
            matches = [abbr for abbr, entry in DATABASES.items() if entry["required"] is (key == "required")]
        elif key in lookup:
            matches = [lookup[key]]
        else:
            unknown.append(name)
            continue
        for abbr in matches:
            for item in [abbr, *DATABASES[abbr].get("with", [])]:
                if item not in selected:
                    selected.append(item)

    if unknown:
        choices = f"Choose from: {', '.join(DATABASES)} (or their file names), or 'required', 'optional', 'taxpurge', 'all'. Run 'AAFTF database' to list them."
        raise ValueError(f"unknown database(s): {', '.join(unknown)}\n{choices}")
    return selected
