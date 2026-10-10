"""AAFTF: Automatic Assembly For The Fungi."""

import logging
from importlib.metadata import PackageNotFoundError, version

__all__ = ["__version__"]


try:
    __version__ = version("AAFTF")
except PackageNotFoundError:
    # imported from a source checkout that was never installed
    __version__ = "0.0.0+unknown"

# Library code only logs; aaftf.main.main() configures handlers via setup_logging().
logging.getLogger(__name__).addHandler(logging.NullHandler())
