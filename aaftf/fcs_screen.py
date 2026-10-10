"""Run NCBI routines to identify contaminant and vectior contigs.

The contaminants are presumably sequences that were not screened out
in the filter step.

This uses NCBI fcs tool for screening which relies on a singularity
engine installed

The default libraries for screening are located in resources.py
and include common Euk, Prok, and MITO contaminants.
"""

import logging
import shutil
from pathlib import Path
from typing import Any

from aaftf.resources import FCSADAPTOR
from aaftf.utility import cleanup_workdir, close_step_log, make_workdir, require_databases, run_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(
    infile: str,
    outfile: str,
    container_engine: str = "singularity",
    workdir: str | None = None,
    image: str | None = None,
    prok: bool = False,
    fcs_script: str | None = None,
    debug: bool = False,
    **kwargs: Any,
) -> None:
    """Screen and trim adaptor/vector sequence with NCBI FCS-adaptor.

    Runs ``run_fcsadaptor.sh`` in a singularity/apptainer or docker container, moves the cleaned
    FASTA to ``outfile``, prints the adaptor report and saves it as
    ``{outfile}.fcs_adaptor_report.txt``.

    Args:
        infile: Assembly FASTA to screen.
        outfile: Output cleaned FASTA.
        container_engine: ``singularity`` or ``docker``.
        workdir: Working directory; a temporary one is created when None.
        image: Singularity image file or docker image reference; defaults to the installed
            image (singularity) or the pinned FCS-adaptor release (docker).
        prok: Screen as prokaryote (``--prok``) instead of eukaryote (``--euk``).
        fcs_script: Path to ``run_fcsadaptor.sh``; found on PATH or in the database when None.
        debug: Show tool output and keep the work directory.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ...); ignored.

    Raises:
        FileNotFoundError: If the chosen container engine is not on PATH.
        ValueError: If ``container_engine`` is not ``singularity`` or ``docker``.
    """
    containerengine = container_engine
    infilename = Path(infile).resolve().name
    tax = "--euk"
    if prok:
        tax = "--prok"
    workdir, custom_workdir = make_workdir(workdir, "fcs_screen")
    try:
        # the wrapper script and singularity image come from `AAFTF database` unless given/found on PATH
        fcsexe = fcs_script or shutil.which("run_fcsadaptor.sh")
        needed = ([] if fcsexe else ["fcs_script"]) + (["fcs_image"] if containerengine == "singularity" and image is None else [])
        found = dict(zip(needed, require_databases(needed, hint="or pass --fcs_script PATH / --image PATH")))
        fcsexe = fcsexe or found["fcs_script"]

        if containerengine == "singularity":
            image = image or found["fcs_image"]
            if shutil.which("singularity") is None and shutil.which("apptainer") is None:
                raise FileNotFoundError("--container_engine singularity requires 'singularity' or 'apptainer' on PATH.")
        elif containerengine == "docker":
            # docker image reference (registry:tag), not a local file; docker itself
            # resolves/pulls it, so no download step is needed here.
            if image is None:
                image = FCSADAPTOR["DOCKERIMAGE"] % (FCSADAPTOR["VERSION"])
            if shutil.which("docker") is None:
                raise FileNotFoundError("--container_engine docker requires 'docker' on PATH.")
        else:
            raise ValueError(f"Unknown --container_engine {containerengine}; use singularity or docker")

        cmd = [fcsexe, "--fasta-input", infile, "--output-dir", workdir, tax, "--container-engine", containerengine, "--image", image]
        result = run_cmd(cmd, debug)

        cleanresult = str(Path(workdir, "cleaned_sequences", infilename))
        if result.returncode != 0 or not Path(cleanresult).is_file():
            raise RuntimeError(f"FCS-adaptor failed (exit {result.returncode}); no {cleanresult} was written (rerun with -v to see its output)")
        if debug:
            logger.info(f"copy from: {cleanresult} -> {outfile}")
        Path(cleanresult).rename(outfile)
        fcsreport = str(Path(workdir, "fcs_adaptor_report.txt"))
        with open(fcsreport) as fh:
            logger.info("FCS report:")
            for line in fh:
                print(line, end="")
        # make a copy of the report to show
        Path(fcsreport).rename(outfile + ".fcs_adaptor_report.txt")
        # cleanup after running
        cleanup_workdir(workdir, debug, custom_workdir)
    finally:
        close_step_log()
