"""Download and cache AAFTF reference databases.

``AAFTF database`` with no arguments lists every database in
``aaftf.resources.DATABASES``: its abbreviation, file, size, which subcommands
use it, and where it is stored. Naming databases (by abbreviation or file name,
or ``all``) downloads them into the database folder (see
``utility.db_dirs``), so later AAFTF commands reuse the cached files.
"""

import logging
import os
import urllib.request
from pathlib import Path
from typing import Any

from aaftf.resources import DATABASES, TAXPURGE_DEFAULT_DATABASES
from aaftf.utility import db_dirs, db_file, db_write_dir, download_file, find_db_file, open_url, warn_if_home_cache

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(databases: list[str] | None = None, force: bool = False, **kwargs: Any) -> None:
    """Execute the ``database`` subcommand.

    Args:
        databases: Database abbreviations or file names (case-insensitive), or ``required``
            (what the default pipeline needs), ``optional`` (the rest), ``taxpurge`` (its default
            sourmash databases) or ``all``. When empty,
            list the databases instead of downloading.
        force: Re-download even if a copy already exists (into the first
            writable database folder).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``debug``, ...); ignored.

    Raises:
        ValueError: If a database name is not recognised.
        RuntimeError: If any download fails.
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
        dest = db_file(entry["filename"], force=force)
        already_present = Path(dest).exists() and not force
        try:
            download_file(entry["url"], dest, force=force)
        except Exception as e:
            errors.append(f"{abbr} ({entry['filename']}): {e}")
            continue
        if entry.get("executable") and not already_present:
            os.chmod(dest, 0o555)

    if errors:
        hint = "Re-run the same command later; files already downloaded are skipped unless you pass --force."
        raise RuntimeError("Some downloads failed:\n" + "\n".join(errors) + f"\n{hint}")
    logger.info("Download complete. Run 'AAFTF database' to see where each database is stored.")


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
            local_path = find_db_file(entry["filename"])
            size: int | None
            if local_path:
                size = Path(local_path).stat().st_size
                downloaded_total += size
                size_label, location = _human_size(size), str(Path(local_path).parent)
            else:
                size = _remote_size(entry["url"])
                if size is None:
                    pending_unknown = True
                else:
                    pending_total += size
                size_label, location = ("unknown" if size is None else _human_size(size)), "not downloaded"
            rows.append((abbr, entry["filename"], size_label, entry["used_by"], location))

    # Any other files present that aren't databases (e.g. leftover reports)
    known_files = {entry["filename"] for entry in DATABASES.values()}
    other_start = len(rows)
    for folder in folders:
        if not folder.is_dir():
            continue
        for f in sorted(f for f in folder.rglob("*") if f.is_file() and f.name not in known_files):
            size = f.stat().st_size
            downloaded_total += size
            rows.append(("-", str(f.relative_to(folder)), _human_size(size), "-", f"{folder} (other file)"))

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
