"""Assess quality statistics of a genome assembly.

This simply gives GC%, N50, L50, Min, Max
contig statistics.
"""

import logging
import re
from pathlib import Path
from typing import Any, TextIO

from Bio import SeqIO
from Bio.Seq import Seq

from aaftf.utility import COMPLEMENT, calc_nx, open_maybe_gz

__all__ = ["run", "genome_asm_stats", "find_telomere", "make_regex_revcomp"]


logger = logging.getLogger(__name__)


def run(input: str, report: str | None = None, telomere_monomer: str = "TAAC{3,5}", telomere_n_repeat: int = 2, telomere_window: int = 200, **kwargs: Any) -> None:
    """Print assembly statistics, including telomere counts, and optionally save them.

    Args:
        input: Assembly FASTA (optionally gzipped).
        report: File to also write the report to.
        telomere_monomer: Telomere monomer regex.
        telomere_n_repeat: Minimum monomer matches at a contig end to call a telomere.
        telomere_window: Number of bp scanned at each contig end.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.

    Raises:
        FileNotFoundError: If ``input`` does not exist.
    """
    if not Path(input).is_file():
        raise FileNotFoundError(f"assembly file not found: {input}")
    if report:
        with open(report, "w") as output_handle:
            genome_asm_stats(input, output_handle, telomere_monomer, telomere_n_repeat, telomere_window)
    else:
        genome_asm_stats(input, None, telomere_monomer, telomere_n_repeat, telomere_window)


def genome_asm_stats(fasta_file: str, output_handle: TextIO | None, telomere_repeat: str, n_minimum: int, telomere_window: int = 200) -> None:
    """Calculate genome assembly statistics, print them and optionally write them to a handle.

    Reports contig count, total length, min/max/median/mean, L50/N50, L90/N90, GC%, N gaps, soft
    masking and telomere counts.

    Args:
        fasta_file: Assembly FASTA (optionally gzipped).
        output_handle: Open handle the report is also written to, or None.
        telomere_repeat: Telomere monomer regex.
        n_minimum: Minimum monomer matches at a contig end to call a telomere.
        telomere_window: Number of bp scanned at each contig end.

    Raises:
        ValueError: If the assembly has no sequences, or only empty ones.
    """
    lengths = []
    gc_bases = 0
    total_ns = 0
    n_gap_count = 0
    total_masked = 0
    telomere_stats = {"TELOMERE FWD": 0, "TELOMERE REV": 0, "T2T SCAFFOLDS": 0}
    with open_maybe_gz(fasta_file) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            lengths.append(len(record))
            forward, reverse = find_telomere(record.seq, telomere_repeat, n_minimum, telomere_window)
            if forward:
                telomere_stats["TELOMERE FWD"] += 1
            if reverse:
                telomere_stats["TELOMERE REV"] += 1
            if forward and reverse:
                telomere_stats["T2T SCAFFOLDS"] += 1

            seq_str = str(record.seq)
            total_masked += sum(1 for c in seq_str if c.islower())
            seq_upper = seq_str.upper()
            gc_bases += sum(seq_upper.count(x) for x in ["G", "C", "S"])
            total_ns += seq_upper.count("N")
            n_gap_count += len(re.findall(r"N+", seq_upper))

    lengths.sort()
    total_len = sum(lengths)
    if total_len == 0:
        raise ValueError(f"{fasta_file} contains no sequence")
    gc = 100.0 * (gc_bases / total_len)
    n50, l50 = calc_nx(lengths, 0.5)
    n90, l90 = calc_nx(lengths, 0.9)
    report = f"Assembly statistics for: {fasta_file}\n"
    report += f"{'CONTIG COUNT':>15}  =  {len(lengths)}\n"
    report += f"{'TOTAL LENGTH':>15}  =  {total_len}\n"
    report += f"{'MIN':>15}  =  {lengths[0]}\n"
    report += f"{'MAX':>15}  =  {lengths[-1]}\n"
    report += f"{'MEDIAN':>15}  =  {lengths[int(len(lengths) / 2)]}\n"
    report += f"{'MEAN':>15}  =  {total_len / len(lengths):.2f}\n"
    report += f"{'L50':>15}  =  {l50}\n"
    report += f"{'N50':>15}  =  {n50}\n"
    report += f"{'L90':>15}  =  {l90}\n"
    report += f"{'N90':>15}  =  {n90}\n"
    report += f"{'GC%':>15}  =  {gc:.2f}\n"
    report += f"{'N GAP COUNT':>15}  =  {n_gap_count}\n"
    report += f"{'TOTAL N BASES':>15}  =  {total_ns}\n"
    if total_masked > 0:
        report += f"{'BASES MASKED':>15}  =  {total_masked}\n"
        report += f"{'PERCENT MASKED':>15}  =  {100.0 * total_masked / total_len:.2f}\n"
    for f in sorted(telomere_stats):
        report += f"{f:>15}  =  {telomere_stats[f]}\n"

    print(report)
    if output_handle:
        output_handle.write(report)


def find_telomere(seq: Seq | str, monomer: str = "TAACCC", min_copies: int = 2, window_size: int = 200) -> tuple[bool, bool]:
    """Check whether each end of a sequence contains telomere repeats.

    Counts non-overlapping matches of the monomer or its reverse complement (case-insensitive)
    in the first and last ``window_size`` bp; sequences no longer than two windows are split in
    half instead. Based on find_telomeres.py from Markus Hiltunen
    (https://github.com/markhilt/genome_analysis_tools).

    Args:
        seq: Nucleotide sequence.
        monomer: Telomere monomer regex.
        min_copies: Minimum matches needed at an end.
        window_size: Number of bp scanned at each end.

    Returns:
        ``(start_has_telomere, end_has_telomere)``.
    """
    fwd_pattern = monomer.strip()
    rev_pattern = make_regex_revcomp(fwd_pattern)

    combined_regex = re.compile(rf"(?:{fwd_pattern}|{rev_pattern})", re.IGNORECASE)

    seq_len = len(seq)

    if seq_len <= window_size * 2:
        mid = seq_len // 2
        start_seq = str(seq[:mid])
        end_seq = str(seq[mid:])
    else:
        start_seq = str(seq[:window_size])
        end_seq = str(seq[-window_size:])

    start_count = len(combined_regex.findall(start_seq))
    end_count = len(combined_regex.findall(end_seq))

    forward = start_count >= min_copies
    reverse = end_count >= min_copies

    return forward, reverse


def make_regex_revcomp(pattern: str) -> str:
    """Reverse complement a DNA sequence or simple regex.

    Handles single bases and ``[...]`` classes with optional ``{m,n}``/``+``/``*``/``?``
    quantifiers, which stay attached to their base; other symbols pass through. Based on
    find_telomeres.py from Markus Hiltunen (https://github.com/markhilt/genome_analysis_tools).

    Args:
        pattern: Sequence or regex to reverse complement.

    Returns:
        The reverse-complemented pattern.
    """
    token_re = re.compile(r"(\[[^\]]+\]|\w)(?:(\{[0-9,]+\}|\+|\*|\?))?")

    tokens = []
    pos = 0
    while pos < len(pattern):
        m = token_re.match(pattern, pos)
        if not m:
            # Pass through any standalone symbol (e.g. delimiters)
            tokens.append((pattern[pos], ""))
            pos += 1
            continue

        base_unit, quantifier = m.group(1), m.group(2) or ""

        # Complement internal characters if it's a bracketed class
        if base_unit.startswith("[") and base_unit.endswith("]"):
            inner = base_unit[1:-1]
            if inner.startswith("^"):
                comp_inner = "^" + inner[1:].translate(COMPLEMENT)
            else:
                comp_inner = inner.translate(COMPLEMENT)
            comp_base = f"[{comp_inner}]"
        else:
            comp_base = base_unit.translate(COMPLEMENT)

        tokens.append((comp_base, quantifier))
        pos = m.end()

    # Invert token order while keeping quantifier bound to its base
    return "".join(base + quant for base, quant in reversed(tokens))
