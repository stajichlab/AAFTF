"""Version reporting for the AAFTF module.

The authoritative version is derived from git release tags (``vX.Y.Z``) by
hatch-vcs at build time, which writes the resolved string into
``AAFTF/_version.py``. Resolution order at runtime:

1. ``git describe`` when running from a source checkout (``.git`` present),
   so the live working tree wins over any stale baked file or installed
   metadata.
2. ``AAFTF/_version.py`` written by hatch-vcs into built/installed packages.
3. Installed package metadata (``pip``/``conda`` installs without a baked file).
4. A static fallback so imports never fail.
"""

from os.path import dirname, isdir, join
from subprocess import DEVNULL, CalledProcessError, check_output

PREFIX = "v"

# Static fallback, used only when no build-time version file, installed
# package metadata, or git checkout is available. Keep roughly in sync with
# the latest release tag; it is not the source of truth.
__fallback_version__ = "0.6.2"


def _version_from_file():
    """Return the version baked in by hatch-vcs, if present."""
    try:
        from AAFTF._version import __version__ as v  # type: ignore

        return v
    except Exception:
        return None


def _version_from_metadata():
    """Return the installed-package version from importlib.metadata, if present."""
    try:
        from importlib.metadata import PackageNotFoundError, version

        return version("AAFTF")
    except PackageNotFoundError:
        return None
    except Exception:
        return None


def _normalize_tag(tag):
    """Strip the ``v`` prefix and normalize a tag to a PEP 440 release string."""
    if tag.startswith(PREFIX):
        tag = tag[len(PREFIX) :]
    try:
        # Normalize e.g. "0.6.1-alpha2" -> "0.6.1a2" when packaging is available.
        from packaging.version import Version

        return str(Version(tag))
    except Exception:
        return tag


def _version_from_git():
    """Build a PEP 440 version from git tags + short hash for a source checkout."""
    d = dirname(dirname(__file__))
    if not isdir(join(d, ".git")):
        return None
    try:
        described = (
            check_output(
                ["git", "describe", "--tags", "--long", "--always", "--dirty", "--match", f"{PREFIX}[0-9]*"],
                cwd=d,
                stderr=DEVNULL,
            )
            .decode()
            .strip()
        )
    except (CalledProcessError, OSError):
        return None

    dirty = described.endswith("-dirty")
    if dirty:
        described = described[: -len("-dirty")]

    # With --long the form is "<tag>-<distance>-g<hash>"; with --always and no
    # tags it is just "<hash>".
    parts = described.rsplit("-", 2)
    if len(parts) == 3 and parts[2].startswith("g"):
        tag, distance, ghash = parts
        base = _normalize_tag(tag)
        short = ghash[1:]  # drop leading 'g'
        if int(distance) == 0 and not dirty:
            return base
        local = short + (".dirty" if dirty else "")
        return f"{base}.dev{distance}+{local}"

    # No tags reachable: just the abbreviated hash (already a node string).
    local = described + (".dirty" if dirty else "")
    return f"0.0.0+{local}"


def get_version():
    """Return the best available version string.

    In a source checkout the live ``git describe`` wins so the working tree's
    state (including ``-dirty``) is always reflected, even if a stale
    ``AAFTF/_version.py`` was left behind by a previous local build. When
    there is no ``.git`` (or no git binary), ``_version_from_git`` returns
    ``None`` and resolution falls through to the baked file / installed
    metadata used by real installs.
    """
    return _version_from_git() or _version_from_file() or _version_from_metadata() or __fallback_version__


__version__ = get_version()


if __name__ == "__main__":
    print(get_version())
