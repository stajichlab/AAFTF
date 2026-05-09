"""Download and cache AAFTF reference databases.

This module provides a one-time download command that downloads all (or a
user-selected subset) of the reference databases listed in
``AAFTF.resources`` into a persistent ``$AAFTF_DB`` folder.  Subsequent
AAFTF commands will reuse these cached files instead of re-downloading
them on first use.
"""

import os
import shutil
import sys
import urllib.request

from AAFTF.resources import FCSADAPTOR, Contaminant_Accessions, DB_Links
from AAFTF.utility import SafeRemove, status


class _Redirect308Handler(urllib.request.HTTPRedirectHandler):
    """Extend urllib's redirect handler to also follow HTTP 308.

    Python < 3.11 does not handle 308 (Permanent Redirect).  The base
    class ``redirect_request()`` hard-codes the allowed set to
    {301,302,303,307} and raises HTTPError for anything else, so we must
    override both that method and add the http_error_308 dispatcher.
    """

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if code == 308:
            code = 307  # method-preserving permanent redirect
        return super().redirect_request(req, fp, code, msg, headers, newurl)

    def http_error_308(self, req, fp, code, msg, headers):
        return self.http_error_302(req, fp, code, msg, headers)


_opener = urllib.request.build_opener(_Redirect308Handler())


def _download(url, dest, force=False):
    """Download ``url`` to ``dest`` if it does not already exist.

    Writes to a temporary file first and renames atomically on success so
    that an interrupted download never leaves a partial file that would be
    mistaken for a complete one on the next run.

    Args:
        url: Remote URL to download.
        dest: Local file path to write.
        force: If True, re-download even if ``dest`` exists.

    Returns:
        The absolute path to the downloaded file.
    """
    if os.path.exists(dest) and not force:
        status(f"  Already present: {dest}")
        return dest

    status(f"  Downloading {os.path.basename(dest)} ...")
    os.makedirs(os.path.dirname(dest), exist_ok=True)

    tmp = dest + ".tmp"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "AAFTF/1.0"})
        with _opener.open(req, timeout=300) as response:
            final_url = response.geturl()
            if final_url != url:
                status(f"  Redirected to {final_url}")
            with open(tmp, "wb") as outfh:
                shutil.copyfileobj(response, outfh)
        os.rename(tmp, dest)
    except Exception as e:
        status(f"  ERROR downloading {url}: {e}")
        if os.path.exists(tmp):
            SafeRemove(tmp)
        raise

    status(f"  Saved {dest}")
    return dest


def _download_db_links(db_dir, keys=None, force=False):
    """Download files referenced in ``DB_Links``.

    Args:
        db_dir: Target directory.
        keys: Iterable of ``DB_Links`` keys to download, or ``None``
            for the subset used by *filter* / *vecscreen*.
        force: If True, overwrite existing files.
    """
    if keys is None:
        # Default set needed by filter / vecscreen
        keys = ("UniVec", "CONTAM_EUKS", "CONTAM_PROKS", "MITO")

    for key in keys:
        if key not in DB_Links:
            status(f"  WARNING: unknown DB_Links key '{key}', skipping")
            continue
        status(f"Downloading {key} ...")
        for url_or_meta in DB_Links[key]:
            if isinstance(url_or_meta, dict):
                url = url_or_meta["url"]
                filename = url_or_meta["filename"]
            else:
                url = url_or_meta
                filename = os.path.basename(url)
            dest = os.path.join(db_dir, filename)
            _download(url, dest, force=force)


def _download_contaminants(db_dir, force=False):
    """Download contaminant accessions (e.g. PhiX).

    Args:
        db_dir: Target directory.
        force: If True, overwrite existing files.
    """
    status("Downloading contaminant accessions ...")
    for name, urls in Contaminant_Accessions.items():
        status(f"  {name} ...")
        for url in urls:
            filename = os.path.basename(url)
            dest = os.path.join(db_dir, filename)
            _download(url, dest, force=force)


def _download_sourmash(db_dir, sourdb_type="gbk", force=False):
    """Download sourmash LCA taxonomy databases.

    Args:
        db_dir: Target directory.
        sourdb_type: One of ``gbk``, ``gtdb``, ``gtdbrep``, or ``all``.
        force: If True, overwrite existing files.
    """
    type_map = {
        "gbk": "sourmash_gbk",
        "gtdb": "sourmash_gtdb",
        "gtdbrep": "sourmash_gtdbrep",
    }

    if sourdb_type == "all":
        indices = list(type_map.values())
    else:
        if sourdb_type not in type_map:
            status(f"  ERROR: unknown sourdb_type '{sourdb_type}'. " f"Choose from {list(type_map.keys()) + ['all']}")
            return
        indices = [type_map[sourdb_type]]

    for idx in indices:
        status(f"Downloading sourmash database ({idx}) ...")
        for entry in DB_Links[idx]:
            dest = os.path.join(db_dir, entry["filename"])
            _download(entry["url"], dest, force=force)


def _download_fcs(db_dir, force=False):
    """Download NCBI FCS-adaptor script and SIF image.

    Args:
        db_dir: Target directory.
        force: If True, overwrite existing files.
    """
    status("Downloading NCBI FCS-adaptor resources ...")

    # Wrapper script
    script_url = FCSADAPTOR["EXEURL"] % FCSADAPTOR["VERSION"]
    script_dest = os.path.join(db_dir, "run_fcsadaptor.sh")
    if os.path.exists(script_dest) and not force:
        status(f"  Already present: {script_dest}")
    else:
        _download(script_url, script_dest, force=force)
        os.chmod(script_dest, 0o555)

    # Singularity image
    image_name = FCSADAPTOR["SIFLOCAL"] % FCSADAPTOR["VERSION"]
    image_dest = os.path.join(db_dir, image_name)
    if os.path.exists(image_dest) and not force:
        status(f"  Already present: {image_dest}")
    else:
        image_url = os.path.join(FCSADAPTOR["SIFURL"], FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"])
        _download(image_url, image_dest, force=force)


def run(parser, args):
    """Execute the ``download`` subcommand.

    Downloads reference databases to the ``AAFTF_DB`` directory so that
    later AAFTF commands do not need to fetch them on-the-fly.
    """
    # Resolve database directory
    db_dir = None
    if args.AAFTF_DB:
        db_dir = args.AAFTF_DB
    elif "AAFTF_DB" in os.environ:
        db_dir = os.environ["AAFTF_DB"]
    else:
        status("ERROR: No database directory specified.\n" "  Set the AAFTF_DB environment variable or pass --AAFTF_DB.\n" "  Example:\n" "    export AAFTF_DB=/path/to/aaftf_db\n" "    AAFTF download")
        sys.exit(1)

    db_dir = os.path.abspath(db_dir)
    os.makedirs(db_dir, exist_ok=True)
    status(f"AAFTF database directory: {db_dir}")

    errors = []

    # Core databases (always downloaded unless --skip-core)
    if not args.skip_core:
        try:
            _download_contaminants(db_dir, force=args.force)
        except Exception as e:
            errors.append(f"Contaminant accessions: {e}")
        try:
            _download_db_links(db_dir, force=args.force)
        except Exception as e:
            errors.append(f"DB_Links: {e}")

    # Sourmash taxonomy databases (optional)
    if not args.skip_sourmash:
        try:
            _download_sourmash(db_dir, sourdb_type=args.sourdb_type, force=args.force)
        except Exception as e:
            errors.append(f"Sourmash DB: {e}")

    # NCBI FCS-adaptor resources (optional)
    if not args.skip_fcs:
        try:
            _download_fcs(db_dir, force=args.force)
        except Exception as e:
            errors.append(f"FCS resources: {e}")

    if errors:
        status("\nSome downloads failed:")
        for err in errors:
            status(f"  - {err}")
        status("\nYou can re-run 'AAFTF download' later; already-downloaded " "files will be skipped unless you pass --force.")
        sys.exit(1)

    status("Setup complete. Future AAFTF runs will use cached files from " f"{db_dir}")
