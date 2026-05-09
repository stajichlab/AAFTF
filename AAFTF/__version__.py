"""Version reporting for the AAFTF module."""

import re
from os.path import dirname, isdir, join
from subprocess import DEVNULL, CalledProcessError, check_output

PREFIX = "v"

tag_re = re.compile(rf"\btag: {PREFIX}([0-9][^,]*)\b")
version_re = re.compile("^Version: (.+)$", re.M)

__version__ = "0.6.2"


def get_version():
    """Return the version if it has been injected into the file by git-archive."""
    version = tag_re.search("$Format:%D$")
    if version:
        return version.group(1)

    d = dirname(dirname(__file__))

    if isdir(join(d, ".git")):
        # Get the version using "git describe".
        cmd = f"git describe --tags --match {PREFIX}[0-9]* --dirty --always"
        try:
            version = check_output(cmd.split(), cwd=d, stderr=DEVNULL).decode().strip()[len(PREFIX) :]
        except CalledProcessError:
            return __version__

        # Get the short commit hash (7 characters)
        try:
            short_hash = check_output(["git", "rev-parse", "--short=7", "HEAD"], cwd=d, stderr=DEVNULL).decode().strip()
        except CalledProcessError:
            short_hash = "unknown"

        # PEP 440 compatibility
        if "-" in version:
            if version.endswith("-dirty"):
                # Include short hash for dirty trees
                version = version.replace("-dirty", f".dirty+{short_hash}")
            else:
                # The git describe already includes a short hash, but let's ensure consistency
                # by replacing it with our explicitly requested 7-char hash
                parts = version.split("-")
                if len(parts) >= 3 and parts[-1].startswith("g"):
                    # Replace the existing git hash with our consistent 7-char hash
                    parts[-1] = f"g{short_hash}"
                    version = "-".join(parts)
                version = ".post".join(version.split("-")[:2]) + f"+{short_hash}"
        else:
            # If no additional commits, still add the short commit hash
            version = f"{version}+{short_hash}"

    else:
        # Extract the version from the PKG-INFO file.
        try:
            with open(join(d, "PKG-INFO")) as f:
                match = version_re.search(f.read())
                if match:
                    version = match.group(1)
                else:
                    version = __version__
        except FileNotFoundError:
            # Fallback to hardcoded version if PKG-INFO is not available
            version = __version__

    return version


if __name__ == "__main__":
    print(get_version())
