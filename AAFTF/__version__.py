"""Version reporting for the AAFTF module.

The authoritative version is derived from git release tags (``vX.Y.Z``) by
hatch-vcs at build time, which writes the resolved string into
``AAFTF/_version.py``. The build uses setuptools-scm (wrapped by hatch-vcs)
with the schemes declared in ``pyproject.toml``; the constants ``VERSION_SCHEME``
and ``LOCAL_SCHEME`` below MUST be kept in sync with that file so a source
checkout reports exactly the string a build of the same commit would bake.

Resolution order at runtime:

1. setuptools-scm computed from git, when running from a source checkout
   (``.git`` present) AND setuptools-scm is importable. This uses the same
   engine and schemes as the build, so the result is byte-identical to a
   freshly built wheel of the same commit, and the live working tree wins
   over any stale baked file.
2. A pure-Python ``git describe`` approximation, used in a source checkout
   when setuptools-scm is not installed. It matches the public version and
   structure but the abbreviated hash length / dirty-date suffix may differ
   from setuptools-scm. Install the ``dev`` extra for exact dev versions.
3. ``AAFTF/_version.py`` written by hatch-vcs into built/installed packages
   (real ``pip``/``conda`` installs have no ``.git``, so steps 1-2 are skipped
   and this baked value — itself produced by setuptools-scm — is used).
4. Installed package metadata (installs without a baked file).
5. A static fallback so imports never fail.
"""

from os.path import dirname, isdir, join
from subprocess import DEVNULL, CalledProcessError, check_output

PREFIX = "v"

# Must match [tool.hatch.version].raw-options in pyproject.toml so the runtime
# (source checkout) version is identical to the one hatch-vcs/setuptools-scm
# bakes at build time.
VERSION_SCHEME = "guess-next-dev"
LOCAL_SCHEME = "node-and-date"

# Static fallback, used only when no setuptools-scm, build-time version file,
# git checkout, or installed package metadata is available. Keep roughly in
# sync with the latest release tag; it is not the source of truth.
__fallback_version__ = "0.6.2"


def _repo_root():
    """Return the directory expected to contain ``.git`` for a source checkout."""
    return dirname(dirname(__file__))


def _version_from_scm():
    """Return setuptools-scm's version for this checkout, identical to a build.

    Returns ``None`` when there is no ``.git`` or setuptools-scm is not
    importable, so installs (which carry a baked ``_version.py`` and no git)
    fall through to that file.
    """
    root = _repo_root()
    if not isdir(join(root, ".git")):
        return None
    try:
        from setuptools_scm import get_version as _scm_get_version
    except Exception:
        return None
    try:
        return _scm_get_version(
            root=root,
            version_scheme=VERSION_SCHEME,
            local_scheme=LOCAL_SCHEME,
        )
    except Exception:
        return None


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


def _guess_next_base(base):
    """Bump ``base`` the way setuptools-scm's ``guess-next-dev`` scheme does.

    Increments the prerelease counter when one is present (``0.6.1a2`` ->
    ``0.6.1a3``), otherwise bumps the last release component (``0.6.1`` ->
    ``0.6.2``). Falls back to the unmodified base if ``packaging`` is absent.
    """
    try:
        from packaging.version import Version

        v = Version(base)
        if v.pre is not None:
            letter, num = v.pre
            return f"{v.base_version}{letter}{num + 1}"
        release = list(v.release)
        release[-1] += 1
        nxt = ".".join(str(x) for x in release)
        return f"{v.epoch}!{nxt}" if v.epoch else nxt
    except Exception:
        return base


def _version_from_git():
    """Approximate setuptools-scm from ``git describe`` (no-dependency fallback).

    Mirrors the ``guess-next-dev`` + ``node-and-date`` schemes used at build
    time; the abbreviated hash length and dirty-date detail may differ from
    setuptools-scm, so this is only used when setuptools-scm is unavailable.
    """
    root = _repo_root()
    if not isdir(join(root, ".git")):
        return None
    try:
        described = (
            check_output(
                ["git", "describe", "--tags", "--long", "--always", "--dirty", "--match", f"{PREFIX}[0-9]*"],
                cwd=root,
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
        if int(distance) == 0 and not dirty:
            return base  # exact tag: setuptools-scm returns the tag verbatim
        nextbase = _guess_next_base(base)
        # node-and-date local: "+g<hash>" when clean, "+g<hash>.d<YYYYMMDD>"
        # (current date) when dirty.
        local = ghash
        if dirty:
            from datetime import date

            local = f"{ghash}.d{date.today():%Y%m%d}"
        return f"{nextbase}.dev{distance}+{local}"

    # No tags reachable: just the abbreviated hash (already a node string).
    local = described + (".dirty" if dirty else "")
    return f"0.0.0+{local}"


def get_version():
    """Return the best available version string.

    In a source checkout the live setuptools-scm value (or its git
    approximation) wins so the working tree's state is always reflected, even
    if a stale ``AAFTF/_version.py`` was left behind by a previous local build.
    When there is no ``.git`` (or no git binary), the git-based steps return
    ``None`` and resolution falls through to the baked file / installed
    metadata used by real installs.
    """
    return _version_from_scm() or _version_from_git() or _version_from_file() or _version_from_metadata() or __fallback_version__


__version__ = get_version()


if __name__ == "__main__":
    print(get_version())
