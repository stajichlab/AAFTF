"""Sort FASTA sequences by length and rename their headers."""

import logging
from typing import Any

from Bio.SeqIO.FastaIO import SimpleFastaParser

from aaftf.utility import write_fasta

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(input: str, out: str, minlen: int = 0, name: str = "scaffold", **kwargs: Any) -> None:
    """Sort contigs longest to shortest, drop short ones, and rename them ``{name}_N``.

    Only the first record for each duplicated header is kept.

    Args:
        input: Input FASTA.
        out: Output FASTA.
        minlen: Minimum sequence length to keep.
        name: Prefix for the new sequence names.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.
    """
    logger.info("Sorting sequences by length longest --> shortest")
    all_seqs = {}
    with open(input) as fasta_in:
        for header, seq in SimpleFastaParser(fasta_in):
            if header not in all_seqs:
                if len(seq) >= minlen:
                    all_seqs[header] = seq
    sorted_seqs = sorted(all_seqs.items(), key=lambda item: len(item[1]), reverse=True)
    with open(out, "w") as fasta_out:
        for i, (header, seq) in enumerate(sorted_seqs):
            write_fasta(fasta_out, f"{name}_{i + 1}", seq)

    logger.info(f"Output written to: {out}")
