# AAFTF Development Guidelines

This document provides guidelines for agentic coding agents working on the AAFTF (Automatic Assembly For The Fungi) codebase.

## Project Overview

AAFTF is a Python-based bioinformatics toolkit for automated genome assembly, cleanup, mitochondrial genome assembly, and polishing focused on fungal genomes. The package uses external bioinformatics tools and provides a unified interface for genome assembly workflows.

## Build and Development Commands

### Installation and Setup
```bash
# pixi (recommended; locked versions in pixi.lock). Environments: default (pipeline tools),
# complete (+ optional tools; feature "optional"), dev (+ pytest, pre-commit, ruff; feature "dev")
pixi install -e dev
pixi run -e dev AAFTF -h
pixi run -e complete install-bowtie2   # optional: build bowtie2 from source (AVX2 issue)

# or conda: environment.yml (pipeline tools, AAFTF from PyPI) /
# environment.dev.yml (all tools, AAFTF installed editable from this checkout)
conda env create -f environment.dev.yml && conda activate aaftf-dev

# the package is built with hatchling + hatch-vcs; Python >=3.10,<3.15
pip install -e .
```

### Code Quality and Linting
```bash
# ruff settings are in ruff.toml (target py310)
# Run ruff for code style, linting, and import sorting (replaces flake8, isort, pyupgrade)
ruff check aaftf/
ruff format aaftf/

# Run ruff with automatic fixes
ruff check --fix aaftf/

# Run pydocstyle for documentation standards
pydocstyle --convention=google aaftf/

# Run codespell for spelling checks
codespell aaftf/ --ignore-words-list=nd,reacher,thist,ths,ure,referenc,wile,nin,pid

# All pre-commit hooks (recommended)
pre-commit run --all-files
```

### Testing
```bash
# Unit tests (no external tools needed)
pixi run -e dev pytest tests/ -m unit -v

# All tests (integration tests need the external tools, e.g. from the pixi environment)
pixi run -e dev pytest tests/ -v
```

## Code Style Guidelines

### Import Organization
- Standard library imports first, then third-party, then local AAFTF imports
- Use explicit imports (avoid `from module import *`)
- Group related imports together
- Import AAFTF modules with absolute imports (`from aaftf.utility import ...`)

```python
# Standard library
import os
import sys
import subprocess

# Third-party
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Local AAFTF imports
from aaftf.utility import check_file, run_cmd
from aaftf.resources import DB_LINKS
```

### Function and Variable Naming
- Use `snake_case` for function and variable names
- Use `UPPER_SNAKE_CASE` for constants
- Function names should be descriptive verbs (`run_assembly`, `check_file`, `calculate_stats`)
- Variable names should be clear and concise (`input_file`, `output_dir`, `cpu_count`)

### Documentation and Docstrings
- Use Google-style docstrings and type hints on all function signatures
- Every module should have a module docstring explaining its purpose
- Every public function should have a docstring with description, arguments, and returns

```python
def calculate_n50(contig_lengths):
    """Calculate N50 statistic from a list of contig lengths.

    Args:
        contig_lengths: List of integers representing contig lengths

    Returns:
        Integer N50 value
    """
```

### Error Handling
- Use try/except blocks for external tool execution
- Provide informative error messages with context
- Report progress with a module-level `logger = logging.getLogger(__name__)` (`logger.info`, `logger.warning`, `logger.error`); `aaftf.main.main()` configures it via `setup_logging()` from `-q/--quiet` and `-v/--verbose` (argparse dest `debug`). Don't add "ERROR:"/"WARNING:" prefixes or leading/trailing whitespace to messages — the formatter adds the level name
- Handle file existence checks with `check_file()` utility

### Code Structure
- Package layout: `aaftf/` (one module per subcommand: `trim.py`, `assemble.py`, ...; `_menu.py` holds every subcommand's argparse parser; `main.py` is the entry point; `pipeline.py`; `utility.py`; `resources.py`; `data/`), `tests/` (pytest, `test_<module>.py`, shared fixtures in `conftest.py`), `docs/` (Sphinx)
- Common utilities are in `utility.py`
- Resource URLs and constants are in `resources.py`
- Each subcommand module has a `run(<CLI dest names>..., **kwargs)` function; its parser is a `<name>_menu()` in `aaftf/_menu.py` that binds it with `set_defaults(func=<module>.run)` (see CLAUDE.md)

### CLI Argument Parsing
- Use `argparse` for command-line interface
- Follow existing patterns for common arguments (`-c/--cpus`, `-m/--memory`, `-o/--out`, `-1/--read1`, `-2/--read2`, `-w/--workdir`); path arguments use `type=str`
- Finish every subcommand's "optional arguments" group with `add_verbosity_args(optional)` (`-q/--quiet`, `-v/--verbose` with dest `debug`); `pipe` is a `run()`-only parameter, not a CLI option
- Provide help text for all arguments
- Use consistent naming across subcommands

### External Tool Integration
- Check tool availability with `require_tools([...])` (raises `FileNotFoundError` naming the missing tools); `AAFTF dependency` reports every tool
- Run external commands with `run_cmd()` / `execute()` (or `align_to_sorted_bam()` for aligner → sorted, indexed BAM)

### File I/O Patterns
- Support both compressed (.gz) and uncompressed files
- Use `check_file()` to validate input files
- Use `safe_remove()` for file cleanup. Create a subcommand's working directory with `workdir, custom_workdir = make_workdir(workdir, "<name>")` and remove it with `cleanup_workdir(workdir, debug, custom_workdir)` (it never deletes a user-supplied `--workdir`)
- Reference databases are listed in `aaftf.resources.DATABASES` and downloaded **only** by `AAFTF database`. Other subcommands get them with `require_databases([...])`, which returns the stored paths or raises `FileNotFoundError` naming the `AAFTF database ...` command to run — never download them from a subcommand. (Only `filter`'s `-a/--screen_accessions` / `-u/--screen_urls` fetch sequences themselves.)
- Every module follows the same top-level layout: (1) module docstring; (2) imports (stdlib, third-party, local; conditional `try: import ...` blocks right after them); (3) `__all__` listing the public constants, classes and functions; (4) constants and module-level variables, public first, then private ones and shared state (`logger`, `_cache`, ...); (5) public ABCs/Protocols; (6) public classes; (7) public functions — `run` first in subcommand modules, the rest in the order they are called; (8) private classes; (9) private helpers, also in call order; (10) the `if __name__ == "__main__":` block. Module-level code must not call anything defined further down (e.g. build objects lazily in a function, as `open_url()` does).
- Other shared helpers in `aaftf/utility.py` — use these instead of re-implementing: `download_file()` (never leaves a partial file behind; for `database` and filter's -a/-u), `db_file()` (database file location across `$AAFTF_DB` folders / `~/.cache/aaftf`), `require_tools()`, `run_cmd()` / `execute()`, `align_to_sorted_bam()`, `concat_files()`, `open_maybe_gz()`, `next_step_name()`, `basename_from_reads()`, `filter_fasta()` / `write_fasta()`, `calc_nx()`
- Never call `sys.exit()` in a subcommand or helper — raise instead, and let `aaftf.main.main()` log the message and return the exit code (the console script passes it to `sys.exit`): `FileNotFoundError` for a missing input file, tool or database (exit 2), `ValueError` for bad or missing arguments, `RuntimeError` when an external tool fails or produces no output (both exit 1); Ctrl-C returns 130. Put the whole message (including any hint) in the exception rather than logging it first, so it is not printed twice. The traceback is shown with `-v/--verbose`.
- Handle file paths with `os.path` operations for cross-platform compatibility

### Logging and Output
- Use `logger.info/warning/error` for user-facing messages
- Use `print_cmd()` to show commands being executed
- Provide progress feedback for long-running operations
- Use debug flags for verbose output during development

### Performance Considerations
- Respect CPU and memory limits specified by user
- Use streaming operations for large files when possible
- Implement efficient file parsing with BioPython iterators
- Consider external tool memory requirements in parameter defaults

### Version Handling
- `aaftf/__init__.py` resolves `AAFTF.__version__` via `importlib.metadata.version("AAFTF")`, falling back to `"0.0.0+unknown"` if the package isn't installed
- The version string itself is set at build time by `hatch-vcs` from git tags (see `[tool.hatch.version]` in `pyproject.toml`)
- Main application (`main.py`) imports it with `from aaftf import __version__`
- Version is displayed both via `--version` flag and at application startup
- Maintains PEP 440 compatibility for Python packaging standards

## Testing Guidelines

Since AAFTF primarily integrates external bioinformatics tools, testing focuses on:
- pytest tests in `tests/`, marked `unit` (no external tools) or `integration` (needs the external tools)
- Mocked subprocess commands and files written (assert on side effects, not on the caller's kwargs)
- CLI argument parsing validation
- File format compatibility (FASTA/FASTQ parsing)

## Common Patterns

### Subcommand Module Structure
```python
"""Module description."""

import logging
from aaftf.utility import check_file, run_cmd

logger = logging.getLogger(__name__)


def run(input, debug=False, pipe=False, **kwargs):
    """Main entry point for subcommand."""
    # Validate inputs
    if not check_file(input):
        raise FileNotFoundError(f"Input file not found: {input}")

    logger.info("Starting processing...")

    # Execute external tools (prints the command; stderr shown only with -v/--verbose)
    run_cmd(["tool", "--input", input], debug)

    logger.info("Processing complete")
```

This project follows bioinformatics best practices with emphasis on reproducibility, proper tool integration, and clear user feedback.
