#!/usr/bin/env python3
"""AAFTF main framework for submodule run."""

# note structure of code taken from poretools
# https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import argparse as ap
import logging
import sys

# AAFTF imports
from aaftf import __version__
from aaftf._menu import register_subcommands
from aaftf.utility import CustomHelpFormatter, available_cpus, finish_logging, get_ram, setup_logging

__all__ = ["main"]


logger = logging.getLogger("aaftf.main")


def main() -> int:
    """Parse the command line, run the chosen subcommand, and return a shell exit code.

    Returns:
        0 on success, 1 on any error, 2 when a required file, tool or database is
        missing (FileNotFoundError), and 130 when interrupted with Ctrl-C.
    """
    #########################################
    # create the top-level parser
    #########################################
    parser = ap.ArgumentParser(prog="AAFTF", usage="%(prog)s [-h] [-V] <command> [options]", formatter_class=CustomHelpFormatter)
    parser.add_argument("-V", "--version", help="Installed AAFTF version", action="version", version=f"%(prog)s {__version__}")

    #########################################
    # create the individual tool parsers
    #########################################
    register_subcommands(parser)

    # process the arguments now
    # if no args then print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args()

    # each subcommand parser binds its own subtool's run() via
    # parser_x.set_defaults(func=<module>.run) in aaftf/_menu.py; a bare
    # "AAFTF" invocation with unrecognized/no subcommand leaves func unset.
    if not hasattr(args, "func"):
        parser.print_help(sys.stderr)
        return 1

    setup_logging(debug=getattr(args, "debug", False), quiet=getattr(args, "quiet", False))
    try:
        logger.info(f"Running AAFTF v{__version__}")
        _cap_to_available(args)
        args.func(**vars(args))
    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 130
    except BrokenPipeError:  # output piped into e.g. `head`, which closed early
        return 0
    except FileNotFoundError as e:
        logger.error("An error occurred: %s", e, exc_info=True)  # traceback shown with -v, always in the log file
        return 2
    except Exception as e:
        logger.error("An error occurred: %s", e, exc_info=True)
        return 1
    finally:
        finish_logging(args.command, getattr(args, "debug", False))
    return 0


def _cap_to_available(args: ap.Namespace) -> None:
    """Lower -c/--cpus and -m/--memory to what this machine/job can provide, with a warning.

    ``args`` is modified in place; options the subcommand lacks (or left unset) are ignored.

    Args:
        args: Parsed command-line arguments.
    """
    cpus = getattr(args, "cpus", None)
    if cpus and cpus > (avail_cpus := available_cpus()):
        logger.warning(f"-c/--cpus {cpus} is more than the {avail_cpus} CPUs available to this job; using {avail_cpus}")
        args.cpus = avail_cpus

    memory = getattr(args, "memory", None)
    if memory is not None and memory > (avail_ram := get_ram()):
        capped = max(int(avail_ram), 1)
        logger.warning(f"-m/--memory {memory} GB is more than the {avail_ram:g} GB of RAM available; using {capped} GB")
        args.memory = capped


if __name__ == "__main__":
    sys.exit(main())
