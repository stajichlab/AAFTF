"""Check that AAFTF's external tool and Python package dependencies are installed."""

import importlib
import logging
import shutil
from pathlib import Path
from typing import Any

__all__ = ["REQUIRED_TOOLS", "OPTIONAL_TOOLS", "REQUIRED_PYTHON_PACKAGES", "OPTIONAL_PYTHON_PACKAGES", "run"]


# external tools needed to run ``AAFTF pipeline_short``, ``pipeline_long`` and ``pipeline_hybrid`` with their default settings (tool -> where it is used)
REQUIRED_TOOLS = {
    "bbduk.sh": "trim, filter",
    "reformat.sh": "trim, filter (paired reads), mito (--subsample)",
    "java": "BBTools (bbduk.sh, reformat.sh), trim (--method trimmomatic)",
    "spades.py": "assemble",
    "blastn": "vecscreen",
    "makeblastdb": "vecscreen",
    "bwa": "covpurge, filter (--aligner bwa), polish, depth (--aligner bwa)",
    "samtools": "covpurge, filter (--aligner bowtie2/bwa/minimap2), polish, depth",
    "minimap2": "rmdup, mito, depth, polish (nextpolish2, racon)",
    "flye": "assemble (--method flye), pipeline_long, pipeline_hybrid",
    "racon": "polish (--method racon), pipeline_long, pipeline_hybrid",
    "polypolish": "polish (--method polypolish), pipeline_hybrid (paired reads)",
    "pypolca": "polish (--method pypolca), pipeline_hybrid (single-end reads)",
    "freebayes": "polish (--method pypolca / polca), pipeline_hybrid (single-end reads)",
    "mosdepth": "depth, pipelines (tag mitochondrial contigs in the final assembly)",
}

# tools only needed by optional steps or non-default options
OPTIONAL_TOOLS = {
    "pigz": "faster read counting (falls back to gzip)",
    "sourmash": "taxpurge (pipelines: only with -p/--phylum)",
    "fastp": "trim (--method fastp)",
    "trimmomatic": "trim (--method trimmomatic)",
    "bowtie2": "filter (--aligner bowtie2)",
    "bowtie2-build": "filter (--aligner bowtie2)",
    "megahit": "assemble (--method megahit)",
    "unicycler": "assemble (--method unicycler), pipeline_hybrid (--method unicycler)",
    "aaftf_polca.sh": "polish (--method polca; installed with AAFTF)",
    "nextPolish2": "polish (--method nextpolish2)",
    "yak": "polish (--method nextpolish2)",
    "NOVOPlasty.pl": "mito",
    "run_fcsadaptor.sh": "fcs_screen (or download it with: AAFTF database fcs_script)",
    "singularity": "fcs_screen (--container_engine singularity)",
    "apptainer": "fcs_screen (--container_engine singularity)",
    "docker": "fcs_screen (--container_engine docker)",
    "run_gx.py": "fcs_gx_purge",
}

# import name -> distribution name
REQUIRED_PYTHON_PACKAGES = {
    "Bio": "biopython",
    "psutil": "psutil",
}

OPTIONAL_PYTHON_PACKAGES = {
    "matplotlib": "matplotlib",  # depth plots
    "sourmash_plugin_branchwater": "sourmash_plugin_branchwater",  # taxpurge (sourmash scripts fastmultigather)
}

logger = logging.getLogger(__name__)


def run(quiet: bool = False, **kwargs: Any) -> None:
    """Check whether AAFTF's external tools and Python packages are installed.

    Tools and packages are checked in two groups: those needed to run the ``pipeline_*`` subcommands with their
    default settings, and those only needed by optional steps (e.g. ``polish``, ``depth``,
    ``mito``) or non-default options (e.g. ``--aligner bowtie2``). An OK/MISSING line is printed
    for each, with its location (or where a missing one is used). Only missing required items
    are an error.

    Args:
        quiet: Print only the OK/MISSING status and name, without locations or where each is used.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``debug``, ...); ignored.

    Raises:
        FileNotFoundError: If a tool or Python package needed by the default pipeline is missing.
    """
    missing_required = _print_tool_table("Checking tools needed by the default pipelines...", _check_tools(REQUIRED_TOOLS), quiet)
    print()
    missing_optional = _print_tool_table("Checking tools for optional steps and non-default options...", _check_tools(OPTIONAL_TOOLS), quiet)
    print()
    missing_packages = _print_package_table("Checking required Python packages...", REQUIRED_PYTHON_PACKAGES)
    missing_optional += _print_package_table("Checking optional Python packages...", OPTIONAL_PYTHON_PACKAGES)
    print()

    _suggest_bowtie2_avx2()
    if missing_optional:
        logger.info(f"NOTE: {len(missing_optional)} optional tool(s)/package(s) missing (only needed for optional steps or non-default options): {', '.join(missing_optional)}")
    if missing_required or missing_packages:
        raise FileNotFoundError(f"missing what the default pipeline needs: {', '.join(missing_required + missing_packages)}")
    logger.info("Everything the default pipeline needs is installed.")


def _suggest_bowtie2_avx2() -> None:
    """Suggest rebuilding bowtie2 for AVX2 when the installed one is the conda build.

    The conda bowtie2 package cannot use AVX2 (it falls back from its x86-64-v3 build at runtime);
    a source build (``install_scripts/pixi_install_bowtie2.sh``) also installs the AVX2
    ``bowtie2-align-s-v256`` binary next to ``bowtie2``, which is how the two are told apart.
    """
    bowtie2 = shutil.which("bowtie2")
    if bowtie2 and not (Path(bowtie2).parent / "bowtie2-align-s-v256").exists():
        logger.info("NOTE: bowtie2 is the conda build, which cannot use AVX2. For AVX2 speed, rebuild it from source: pixi run install-bowtie2 (pixi), or bash install_scripts/pixi_install_bowtie2.sh (conda)")


def _print_tool_table(title: str, results: list[tuple[str, str | None, str]], quiet: bool = False) -> list[str]:
    """Log ``title`` and print an OK/MISSING line for each checked tool.

    Args:
        title: Heading logged before the table.
        results: ``(tool, path, used_by)`` tuples from ``_check_tools``; ``path`` is None if not found.
        quiet: Print only the status and tool name, without its location or the subcommands using it.

    Returns:
        Names of the tools that were not found.
    """
    logger.info(title)
    missing = []
    for tool, path, used_by in results:
        if path:
            print(f"  [OK]      {tool}" if quiet else f"  [OK]      {tool:<20} {path}")
        else:
            print(f"  [MISSING] {tool}" if quiet else f"  [MISSING] {tool:<20} used by: {used_by}")
            missing.append(tool)
    return missing


def _print_package_table(title: str, packages: dict[str, str]) -> list[str]:
    """Log ``title`` and print an OK/MISSING line for each Python package.

    Args:
        title: Heading logged before the table.
        packages: Mapping of import name to distribution name.

    Returns:
        Distribution names of the packages that could not be imported.
    """
    logger.info(title)
    missing = []
    for dist_name, found in _check_python_packages(packages):
        print(f"  [OK]      {dist_name}" if found else f"  [MISSING] {dist_name}")
        if not found:
            missing.append(dist_name)
    return missing


def _check_tools(tools: dict[str, str]) -> list[tuple[str, str | None, str]]:
    """Look up each tool on ``PATH``.

    Args:
        tools: Mapping of executable name to the subcommands that use it.

    Returns:
        ``(tool, path, used_by)`` tuples in the order ``tools`` lists them; ``path`` is None if not on ``PATH``.
    """
    results = []
    for tool, used_by in tools.items():
        path = shutil.which(tool)
        results.append((tool, path, used_by))
    return results


def _check_python_packages(packages: dict[str, str]) -> list[tuple[str, bool]]:
    """Try importing each Python package.

    Args:
        packages: Mapping of import name to distribution name.

    Returns:
        ``(distribution name, importable)`` tuples in the order ``packages`` lists them.
    """
    results = []
    for module_name, dist_name in packages.items():
        try:
            importlib.import_module(module_name)
            results.append((dist_name, True))
        except ImportError:
            results.append((dist_name, False))
    return results
