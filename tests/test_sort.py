"""Unit tests for AAFTF/sort.py.

sort.run() is pure Python (BioPython only) — no external tools required.
"""

import os
from argparse import Namespace

import pytest
from Bio.SeqIO.FastaIO import SimpleFastaParser

from AAFTF.sort import run

from tests.conftest import SEQ1, SEQ2, SEQ3, SMALL_FASTA

pytestmark = pytest.mark.unit


def _read_fasta(path):
    """Return list of (header, seq) from a FASTA file."""
    records = []
    with open(path) as fh:
        for header, seq in SimpleFastaParser(fh):
            records.append((header, seq))
    return records


def _make_args(tmp_path, fasta_file, minlen=0, name="scaffold"):
    return Namespace(
        input=str(fasta_file),
        out=str(tmp_path / "sorted.fasta"),
        minlen=minlen,
        name=name,
        debug=False,
        pipe=True,
    )


# ---------------------------------------------------------------------------
# Basic sorting
# ---------------------------------------------------------------------------

class TestSortOrder:
    def test_output_file_created(self, sort_args):
        run(None, sort_args)
        assert os.path.exists(sort_args.out)

    def test_sorted_longest_first(self, sort_args):
        run(None, sort_args)
        records = _read_fasta(sort_args.out)
        lengths = [len(seq) for _, seq in records]
        # SEQ3=150, SEQ2=80, SEQ1=20 → descending
        assert lengths == sorted(lengths, reverse=True)

    def test_first_record_is_longest(self, sort_args):
        run(None, sort_args)
        records = _read_fasta(sort_args.out)
        assert len(records[0][1]) == 150   # SEQ3

    def test_last_record_is_shortest(self, sort_args):
        run(None, sort_args)
        records = _read_fasta(sort_args.out)
        assert len(records[-1][1]) == 20   # SEQ1

    def test_all_sequences_present(self, sort_args):
        run(None, sort_args)
        records = _read_fasta(sort_args.out)
        assert len(records) == 3


# ---------------------------------------------------------------------------
# Header renaming
# ---------------------------------------------------------------------------

class TestSortRenaming:
    def test_default_prefix(self, sort_args):
        run(None, sort_args)
        records = _read_fasta(sort_args.out)
        headers = [h for h, _ in records]
        assert headers[0] == "scaffold_1"
        assert headers[1] == "scaffold_2"
        assert headers[2] == "scaffold_3"

    def test_custom_prefix(self, tmp_path, fasta_file):
        args = _make_args(tmp_path, fasta_file, name="contig")
        run(None, args)
        records = _read_fasta(args.out)
        headers = [h for h, _ in records]
        assert headers[0] == "contig_1"
        assert headers[-1] == f"contig_{len(records)}"


# ---------------------------------------------------------------------------
# minlen filtering
# ---------------------------------------------------------------------------

class TestSortMinlen:
    def test_minlen_removes_short_contigs(self, tmp_path, fasta_file):
        # SEQ1 is 20 bp; minlen=50 removes it
        args = _make_args(tmp_path, fasta_file, minlen=50)
        run(None, args)
        records = _read_fasta(args.out)
        assert len(records) == 2
        assert all(len(seq) >= 50 for _, seq in records)

    def test_minlen_zero_keeps_all(self, tmp_path, fasta_file):
        args = _make_args(tmp_path, fasta_file, minlen=0)
        run(None, args)
        records = _read_fasta(args.out)
        assert len(records) == 3

    def test_minlen_removes_all_below_threshold(self, tmp_path, fasta_file):
        # All seqs are ≤ 150 bp; minlen=200 removes everything
        args = _make_args(tmp_path, fasta_file, minlen=200)
        run(None, args)
        records = _read_fasta(args.out)
        assert len(records) == 0

    def test_minlen_exact_boundary_included(self, tmp_path, fasta_file):
        # SEQ2 is exactly 80 bp; minlen=80 should keep it
        args = _make_args(tmp_path, fasta_file, minlen=80)
        run(None, args)
        records = _read_fasta(args.out)
        assert any(len(seq) == 80 for _, seq in records)


# ---------------------------------------------------------------------------
# Sequence content preserved
# ---------------------------------------------------------------------------

class TestSortContent:
    def test_sequences_are_not_modified(self, sort_args):
        run(None, sort_args)
        records = _read_fasta(sort_args.out)
        # Collect original sequences by length for comparison
        originals = {len(SEQ3): SEQ3.upper(),
                     len(SEQ2): SEQ2.upper(),
                     len(SEQ1): SEQ1.upper()}
        for header, seq in records:
            expected = originals[len(seq)]
            # softwrap may add newlines within the stored seq string but
            # SimpleFastaParser strips them; just compare content
            assert seq.upper().replace("\n", "") == expected
