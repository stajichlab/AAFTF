"""Run the fcs-gx tool to look for contaminants."""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

from aaftf.utility import check_file, cleanup_workdir, close_step_log, fasta_stats, filter_fasta, make_workdir, next_step_name, print_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)


# logging - we may need to think about whether this has
# separate name for the different runfolder
# flake8: noqa: C901
def run(
    input: str,
    outfile: str,
    workdir: str | None = None,
    taxid: int = 4890,
    db: str = "/my_tmpfs/gxdb/all",
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run NCBI fcs_gx routines to detect and remove contaminant contigs.

    Sequences the FCS-GX report marks ``EXCLUDE`` are dropped; other actions (``TRIM``, ``FIX``,
    ``REVIEW``, ...) are only logged as warnings, since removing part of a contig is not supported
    here. The report is copied next to ``outfile`` as ``<outfile stem>.fcs_gx-taxonomy.tsv``.

    Args:
        input: Assembly FASTA to screen.
        outfile: Output FASTA of retained contigs.
        workdir: Working directory; a temporary one is created when None.
        taxid: NCBI taxonomy ID of the expected organism.
        db: FCS-GX database path prefix (``{db}.gxi`` must exist).
        debug: Print dropped contigs and keep the work directory.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.

    Raises:
        FileNotFoundError: If ``{db}.gxi`` does not exist.
        RuntimeError: If ``run_gx.py`` fails or writes no report.
    """
    workdir, custom_workdir = make_workdir(workdir, "fcs_gx_purge")
    try:
        # parse database locations
        if not db or not Path(f"{db}.gxi").is_file():
            raise FileNotFoundError(f"{db}.gxi not found, needs to have setup the fcs_gx db - https://github.com/ncbi/fcs/wiki/FCS-GX")

        num_seqs, assembly_size = fasta_stats(str(Path(input)))
        logger.info(f"Assembly is {num_seqs:,} contigs and {assembly_size:,} bp")

        # now filter for taxonomy with sourmash lca classify
        logger.info("Running fcs_gx to get taxonomy classification for each contig")

        # python scripts/run_gx.py --bin-dir dist --gx-db /sw/db/gxdb --tax-id 4842 --fasta

        fcsgx_compute = [
            "run_gx.py",
            "--fasta",
            input,
            "--tax-id",
            f"{taxid}",  # taxid is a numeric
            "--gx-db",
            db,
            "--out-dir",
            workdir,
        ]
        print_cmd(fcsgx_compute)
        fcs_log = "fcs_gx.log"
        with open(str(Path(workdir, fcs_log)), "w") as logfile:
            result = subprocess.run(fcsgx_compute, stderr=logfile)

        fname = Path(input).stem
        # output tsv:
        # seq_id	start_pos	end_pos	seq_len	action	div	agg_cont_cov	top_tax_name

        seqs_to_drop = {}

        fcsgx_tsv = str(Path(workdir, f"{fname}.{taxid}.fcs_gx_report.txt"))
        if result.returncode != 0 or not Path(fcsgx_tsv).is_file():
            raise RuntimeError(f"run_gx.py failed (exit {result.returncode}) or did not write {fcsgx_tsv}; see {Path(workdir, fcs_log)}")
        with open(fcsgx_tsv) as fcsgx_out:
            for line in fcsgx_out:
                cols = line.rstrip("\n").split("\t")
                if line.startswith("#") or len(cols) < 5:
                    continue
                action = cols[4].removeprefix("ACTION_")
                if action == "EXCLUDE":
                    seqs_to_drop[cols[0]] = 1
                else:
                    logger.warning(f"FCS-GX reports {action} for {cols[0]}:{cols[1]}..{cols[2]}; not applied, review it manually")

        # drop contigs from taxonomy before calculating coverage
        logger.info(f"Dropping {len(seqs_to_drop)} contigs from fcs-gx taxonomy screen")
        num_seqs, assembly_size = filter_fasta(input, outfile, lambda seq_id: seq_id not in seqs_to_drop)

        if debug:
            print("Contigs dropped due to taxonomy: {:}".format(",".join(seqs_to_drop)))

        logger.info(f"fcs-gx assembly is {num_seqs:,} contigs and {assembly_size:,} bp")
        next_out = next_step_name(outfile, ".rmdup.fasta")

        if check_file(fcsgx_tsv):
            outbase = Path(outfile).name
            basedir = str(Path(outfile).resolve().parent)
            if "." in outbase:
                outbase = outbase.rsplit(".", 1)[0]
            shutil.copy(fcsgx_tsv, str(Path(basedir, f"{outbase}.fcs_gx-taxonomy.tsv")))

        cleanup_workdir(workdir, debug, custom_workdir)

        if not pipe:
            logger.info(f"Your next command might be:\nAAFTF rmdup -i {outfile} -o {next_out}")
    finally:
        close_step_log()
