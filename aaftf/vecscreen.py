"""Trim and split adaptor/vector sequence out of an assembly, the way NCBI FCS-adaptor does.

Contigs are searched with blastn using the VecScreen preset (``VECSCREEN_TASK_ARGS``) against the ``contamdb`` BLAST database, built by
``AAFTF database contamdb`` (and shared with ``filter``) from UniVec_Core (without its
fungal/yeast-derived entries), PhiX and the BBTools lambda and adapter sequences. Each hit is then acted on with FCS-adaptor's rules:

* Terminal hits (within ``TERMINAL_WINDOW`` bp of a contig end) with a raw score >=
  ``TERMINAL_MIN_SCORE`` (or >= ``TIP_MIN_SCORE`` when the hit reaches the very end) are trimmed
  off, together with everything between the hit and that end (ACTION_TRIM).
* Internal hits must be Strong (raw score >= ``INTERNAL_MIN_SCORE``) to be cut out, splitting the
  contig, so weak/moderate matches inside conserved genes never fragment chromosomes.
* A contig with more than ``MAX_VECTOR_FRACTION`` vector sequence, and any piece shorter than
  ``MIN_CONTIG_LEN`` bp, is dropped (ACTION_EXCLUDE).

VecScreen score categories (https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/):
Strong: terminal >= 24, internal >= 30; Moderate: terminal 19-23, internal 25-29;
Weak: terminal 16-18, internal 23-24.
"""

import csv
import logging
from pathlib import Path
from typing import Any

# biopython needed
from Bio import SeqIO

from aaftf.database import require_contamdb
from aaftf.utility import cleanup_workdir, close_step_log, make_workdir, next_step_name, require_tools, run_cmd, write_fasta

__all__ = ["INTERNAL_MIN_SCORE", "MAX_VECTOR_FRACTION", "MIN_CONTIG_LEN", "TERMINAL_MIN_SCORE", "TERMINAL_WINDOW", "TIP_MIN_SCORE", "VECSCREEN_TASK_ARGS", "run"]

# hits starting/ending within this many bp of a contig end are terminal
TERMINAL_WINDOW = 50
# minimum raw score of a terminal hit, and of one reaching the very first/last base
TERMINAL_MIN_SCORE = 19
TIP_MIN_SCORE = 16
# minimum raw score (Strong) of an internal hit
INTERNAL_MIN_SCORE = 30
# contigs (or pieces) shorter than this after trimming, or more than this fraction vector, are dropped
MIN_CONTIG_LEN = 200
MAX_VECTOR_FRACTION = 0.5

_OUTFMT = "6 qaccver saccver pident length qstart qend evalue bitscore score qlen"

# The BLAST library's VecScreen preset (CBlastVecScreenOptionsHandle: reward 1, penalty -5, gap
# open/extend 3/3, DUST on, e-value 700 over NCBI's fixed 1.75e12 search space). The blastn command
# line has no ``-task vecscreen``, so the preset's values are given explicitly; they are what make
# the raw scores comparable with the VecScreen/FCS-adaptor thresholds above.
VECSCREEN_TASK_ARGS = ["-task", "blastn", "-reward", "1", "-penalty", "-5", "-gapopen", "3", "-gapextend", "3", "-dust", "yes", "-evalue", "700", "-searchsp", "1750000000000"]

logger = logging.getLogger(__name__)


def run(
    infile: str,
    outfile: str,
    workdir: str | None = None,
    cpus: int = 1,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Screen an assembly for adaptor/vector sequence and write the trimmed/split assembly.

    Args:
        infile: Input assembly FASTA.
        outfile: Output cleaned FASTA.
        workdir: Working directory; a temporary one is created if None.
        cpus: Number of BLAST threads.
        debug: Keep the working directory when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``); ignored.
    """
    require_tools(["blastn"])
    blastdb = require_contamdb(blast=True)
    workdir, custom_workdir = make_workdir(workdir, "vecscreen")
    try:
        report = str(Path(workdir, Path(outfile).name + ".vecscreen.tab"))
        cmd = ["blastn", *VECSCREEN_TASK_ARGS, "-db", blastdb, "-query", infile, "-outfmt", _OUTFMT, "-num_threads", str(cpus), "-out", report]
        if run_cmd(cmd, debug).returncode != 0:
            raise RuntimeError(f"blastn failed screening {infile} against {blastdb}")

        hits = _read_vector_hits(report)
        n_in, n_out = _write_cleaned(infile, hits, outfile)
        logger.info(f"{n_in:,} contigs screened, {len(hits):,} with adaptor/vector hits; wrote {n_out:,} contigs to {outfile}")

        if not pipe:
            next_out = next_step_name(outfile, ".taxpurge.fasta")
            logger.info(f"Your next command might be:\nAAFTF taxpurge -i {outfile} -o {next_out} --phylum Ascomycota")

        cleanup_workdir(workdir, debug, custom_workdir)
    finally:
        close_step_log()


def _read_vector_hits(report: str) -> dict[str, list[tuple[int, int]]]:
    """Return the hits of a ``_OUTFMT`` BLASTN report that qualify, as 1-based inclusive spans per contig.

    A hit is terminal when it starts within ``TERMINAL_WINDOW`` bp of the 5' end or ends within
    that of the 3' end; it qualifies with a raw score >= ``TERMINAL_MIN_SCORE``, or >=
    ``TIP_MIN_SCORE`` when it covers the first or last base. An internal hit qualifies with a raw
    score >= ``INTERNAL_MIN_SCORE``.

    Args:
        report: BLASTN tab report.

    Returns:
        Mapping of contig id to its qualifying ``(start, end)`` spans.
    """
    hits: dict[str, list[tuple[int, int]]] = {}
    with open(report) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            qaccver, score, qlen = row[0], int(row[8]), int(row[9])
            start, end = sorted((int(row[4]), int(row[5])))
            if start == 1 or end == qlen:
                min_score = TIP_MIN_SCORE
            elif start <= TERMINAL_WINDOW or qlen - end < TERMINAL_WINDOW:
                min_score = TERMINAL_MIN_SCORE
            else:
                min_score = INTERNAL_MIN_SCORE
            if score >= min_score:
                hits.setdefault(qaccver, []).append((start, end))
    return hits


def _vector_regions(spans: list[tuple[int, int]], length: int) -> list[tuple[int, int]]:
    """Merge the vector spans of one contig into the 0-based, half-open regions to remove.

    Overlapping spans, and spans separated by fewer than ``TERMINAL_WINDOW`` bp (VecScreen's
    "segment of suspect origin"), are merged; a region within ``TERMINAL_WINDOW`` bp of an end is
    extended to that end.

    Args:
        spans: 1-based inclusive ``(start, end)`` hit spans.
        length: Contig length.

    Returns:
        Sorted, non-overlapping ``(start, end)`` regions.
    """
    regions: list[list[int]] = []
    for start, end in sorted(spans):
        start -= 1
        if regions and start - regions[-1][1] < TERMINAL_WINDOW:
            regions[-1][1] = max(regions[-1][1], end)
        else:
            regions.append([start, end])
    if regions[0][0] < TERMINAL_WINDOW:
        regions[0][0] = 0
    if length - regions[-1][1] < TERMINAL_WINDOW:
        regions[-1][1] = length
    return [(start, end) for start, end in regions]


def _write_cleaned(infile: str, hits: dict[str, list[tuple[int, int]]], outfile: str) -> tuple[int, int]:
    """Write ``outfile``: ``infile`` with its vector regions trimmed off or cut out.

    A contig keeps its name when trimmed only at its ends; when cut internally, its pieces are
    named ``split<N>_<id>``. Contigs more than ``MAX_VECTOR_FRACTION`` vector, and pieces (or
    contigs) shorter than ``MIN_CONTIG_LEN`` bp, are left out.

    Args:
        infile: Input assembly FASTA.
        hits: Qualifying vector spans per contig id, from ``_read_vector_hits``.
        outfile: Output FASTA.

    Returns:
        Number of contigs read and of sequences written.
    """
    n_in = n_out = 0
    with open(infile) as fh, open(outfile, "w") as out:
        for record in SeqIO.parse(fh, "fasta"):
            n_in += 1
            seq = str(record.seq)
            pieces = [(0, len(seq))]
            if record.id in hits:
                regions = _vector_regions(hits[record.id], len(seq))
                vector_bp = sum(end - start for start, end in regions)
                if vector_bp > MAX_VECTOR_FRACTION * len(seq):
                    logger.info(f"Excluding {record.id}: {vector_bp:,} of {len(seq):,} bp are adaptor/vector")
                    continue
                bounds = [0] + [pos for region in regions for pos in region] + [len(seq)]
                pieces = [(bounds[i], bounds[i + 1]) for i in range(0, len(bounds), 2) if bounds[i + 1] > bounds[i]]
                logger.info(f"{record.id}: removing adaptor/vector regions {regions}, keeping {pieces}")
            for num, (start, end) in enumerate(pieces, start=1):
                if end - start < MIN_CONTIG_LEN:
                    continue
                name = record.id if len(pieces) == 1 else f"split{num}_{record.id}"
                write_fasta(out, name, seq[start:end])
                n_out += 1
    return n_in, n_out
