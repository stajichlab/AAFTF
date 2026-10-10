"""RMDup module step removes smaller contigs redundant with larger ones.

This uses minimap2 to map each small contig against the longer contigs in the
assembly and removes those which are redundant.
"""

import logging
import os
import sys
from pathlib import Path
from typing import Any

from Bio.SeqIO.FastaIO import SimpleFastaParser

from aaftf.utility import calc_nx, cleanup_workdir, close_step_log, filter_fasta, make_workdir, next_step_name, paf_hits, write_fasta

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(
    input: str,
    out: str,
    workdir: str | None = None,
    cpus: int = 1,
    percent_id: float = 95,
    percent_cov: float = 95,
    minlen: int = 500,
    exhaustive: bool = False,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run routines to identify and remove duplicate contigs.

    Contigs are keyed by ID (the first word of the header, as minimap2 reports
    them). Each contig shorter than N75 (or every contig with ``exhaustive``)
    is aligned against all longer contigs; it is dropped if a hit exceeds both
    ``percent_id`` identity and ``percent_cov`` coverage. Contigs shorter than
    ``minlen`` are always dropped.

    Args:
        input: Assembly FASTA.
        out: Output FASTA of retained contigs.
        workdir: Working directory; a temporary one is created when None.
        cpus: Number of minimap2 threads.
        percent_id: Identity threshold (percent) a hit must exceed.
        percent_cov: Query coverage threshold (percent) a hit must exceed.
        minlen: Contigs shorter than this are dropped.
        exhaustive: Check every contig, not only those shorter than N75.
        debug: Log per-contig progress and keep the work directory.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.
    """
    workdir, custom_workdir = make_workdir(workdir, "rmdup")
    try:
        if debug:
            logger.info(f"input={input} out={out} workdir={workdir} cpus={cpus} percent_id={percent_id} percent_cov={percent_cov} minlen={minlen} exhaustive={exhaustive}")
        logger.info("Looping through assembly shortest --> longest searching for duplicated contigs using minimap2")

        # read the assembly once; later lookups and the per-contig query/reference files use it
        lengths = []
        seqs: dict[str, str] = {}
        with open(input) as infile:
            for header, seq in SimpleFastaParser(infile):
                lengths.append(len(seq))
                seqs.setdefault(header.split(None, 1)[0], seq)
        n50, _ = calc_nx(lengths, 0.5)
        n75, _ = calc_nx(lengths, 0.75)
        logger.info(f"Assembly is {len(lengths):,} contigs; {sum(lengths):,} bp; N50 is {n50:,} bp; N75 is {n75:,} bp")

        # (id, length) sorted shortest --> longest
        by_length = sorted(((seq_id, len(seq)) for seq_id, seq in seqs.items()), key=lambda item: item[1])
        if exhaustive:
            n75 = by_length[-1][1]
        to_check = [item for item in by_length if item[1] < n75]
        logger.info(f"Will check {len(to_check):,} contigs for duplication --> those that are < {n75:,} && > {minlen:,}")

        ignore = set()
        prefix = str(os.getpid())
        for i, (seq_id, length) in enumerate(by_length):
            if length < minlen:
                ignore.add(seq_id)
                continue
            if length > n75:
                sys.stdout.write("\n")
                sys.stdout.flush()
                break
            if debug:
                logger.info(f"Working on {seq_id} len={length} remove_tally={len(ignore)}")
            else:
                sys.stdout.write(f"\rProgress: {i} of {len(to_check)}; remove tally={len(ignore):,}; current={seq_id}; length={length}     ")
                sys.stdout.flush()
            longer = [other_id for other_id, _ in by_length[i + 1 :]]
            query, reference = _write_query_and_reference(seqs, seq_id, longer, workdir, prefix)
            if _is_duplicate(query, reference, seq_id, cpus, percent_id, percent_cov, debug):
                ignore.add(seq_id)

        num_seqs, assembly_size = filter_fasta(input, out, lambda seq_id: seq_id not in ignore)
        logger.info(f"Cleaned assembly is {num_seqs:,} contigs and {assembly_size:,} bp")
        next_out = next_step_name(out, ".final.fasta")

        if not pipe:
            logger.info(f"Your next command might be:\nAAFTF sort -i {out} -o {next_out}")

        cleanup_workdir(workdir, debug, custom_workdir)
    finally:
        close_step_log()


def _write_query_and_reference(seqs: dict[str, str], query_id: str, reference_ids: list[str], workdir: str, prefix: str) -> tuple[str, str]:
    """Write ``query_id`` and ``reference_ids`` (from ``seqs``) to separate FASTA files; return their paths.

    Args:
        seqs: Sequences keyed by contig ID.
        query_id: ID of the query contig.
        reference_ids: IDs of the reference contigs.
        workdir: Directory for the files.
        prefix: Filename prefix (``{prefix}query.fasta``, ``{prefix}reference.fasta``).

    Returns:
        ``(query_path, reference_path)``.
    """
    query = str(Path(workdir, f"{prefix}query.fasta"))
    reference = str(Path(workdir, f"{prefix}reference.fasta"))
    with open(query, "w") as qout:
        write_fasta(qout, query_id, seqs[query_id])
    with open(reference, "w") as rout:
        for ref_id in reference_ids:
            write_fasta(rout, ref_id, seqs[ref_id])
    return query, reference


def _is_duplicate(query: str, reference: str, name: str, cpus: int, percent_id: float, percent_cov: float, debug: bool) -> bool:
    """Return True if minimap2 aligns ``query`` to ``reference`` above both identity and coverage thresholds.

    Args:
        query: Query FASTA path.
        reference: Reference FASTA path.
        name: Query contig ID (for logging).
        cpus: Number of minimap2 threads.
        percent_id: Identity threshold (percent) a hit must exceed.
        percent_cov: Query coverage threshold (percent) a hit must exceed.
        debug: Passed to ``paf_hits``.

    Returns:
        Whether any hit exceeds both thresholds.
    """
    cmd = ["minimap2", "-t", str(cpus), "-x", "asm5", "-N5", reference, query]
    for hit in paf_hits(cmd, debug=debug, quiet=True):
        pident = hit.matches / hit.aln_len * 100
        cov = hit.aln_len / hit.query_len * 100
        logger.debug(f"query={hit.query} hit={hit.target} pident={pident:.2f} coverage={cov:.2f}")
        if pident > percent_id and cov > percent_cov:
            logger.debug(f"{name} duplicated: {pident:.0f}% identity over {cov:.0f}% of the contig. length={hit.query_len}")
            return True
    return False
