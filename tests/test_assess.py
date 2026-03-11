"""Unit tests for AAFTF/assess.py.

Tests cover the pure-Python helper functions (revcomp, findTelomere)
and the full genome_asm_stats / run() pipeline against known small
FASTA files.  No external tools are required.
"""

import io
import os

import pytest

from AAFTF.assess import findTelomere, genome_asm_stats, revcomp, run

from tests.conftest import (
    SMALL_FASTA,
    TELO_FWD_ONLY,
    TELO_NONE,
    TELO_REV_ONLY,
    TELO_T2T,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# revcomp  (assess.py version — regex-aware)
# ---------------------------------------------------------------------------

class TestRevcomp:
    """assess.revcomp handles regex bracket metacharacters in the monomer."""

    def test_simple_dna(self):
        assert revcomp("ATCG") == "CGAT"

    def test_all_complement(self):
        assert revcomp("AAAA") == "TTTT"
        assert revcomp("CCCC") == "GGGG"

    def test_regex_monomer(self):
        # revcomp of "TAA[C]+" should be "[G]+TTA"
        rc = revcomp("TAA[C]+")
        assert rc == "[G]+TTA"

    def test_symmetric(self):
        # revcomp(revcomp(seq)) should equal seq for plain DNA
        seq = "ATCGATCG"
        assert revcomp(revcomp(seq)) == seq


# ---------------------------------------------------------------------------
# findTelomere
# ---------------------------------------------------------------------------

class TestFindTelomere:
    MONOMER = "TAA[C]+"
    N = 2

    def _check(self, seq_str):
        from Bio.Seq import Seq
        return findTelomere(Seq(seq_str), self.MONOMER, self.N)

    def test_t2t_both_ends(self):
        fwd, rev = self._check(TELO_T2T)
        assert fwd is True
        assert rev is True

    def test_forward_only(self):
        fwd, rev = self._check(TELO_FWD_ONLY)
        assert fwd is True
        assert rev is False

    def test_reverse_only(self):
        fwd, rev = self._check(TELO_REV_ONLY)
        assert fwd is False
        assert rev is True

    def test_no_telomere(self):
        fwd, rev = self._check(TELO_NONE)
        assert fwd is False
        assert rev is False


# ---------------------------------------------------------------------------
# genome_asm_stats
# ---------------------------------------------------------------------------

class TestGenomeAsmStats:
    MONOMER = "TAA[C]+"
    N = 2

    def test_contig_count(self, fasta_file, capsys):
        genome_asm_stats(str(fasta_file), None, self.MONOMER, self.N)
        out = capsys.readouterr().out
        assert "CONTIG COUNT" in out
        assert "3" in out

    def test_total_length(self, fasta_file, capsys):
        # SEQ1=20 + SEQ2=80 + SEQ3=150 = 250
        genome_asm_stats(str(fasta_file), None, self.MONOMER, self.N)
        out = capsys.readouterr().out
        assert "250" in out

    def test_gc_percent(self, fasta_file, capsys):
        # GC = (10 + 80 + 0) / 250 * 100 = 36.00
        genome_asm_stats(str(fasta_file), None, self.MONOMER, self.N)
        out = capsys.readouterr().out
        assert "36.00" in out

    def test_n50(self, fasta_file, capsys):
        # N50 = 150 (first contig from largest reaches 50 % cumulative)
        genome_asm_stats(str(fasta_file), None, self.MONOMER, self.N)
        out = capsys.readouterr().out
        assert "N50" in out
        # N50 value appears in the report
        assert "150" in out

    def test_l50(self, fasta_file, capsys):
        # L50 = 1
        genome_asm_stats(str(fasta_file), None, self.MONOMER, self.N)
        out = capsys.readouterr().out
        assert "L50" in out
        # L50 value = 1
        assert "1" in out

    def test_writes_to_output_handle(self, fasta_file):
        buf = io.StringIO()
        genome_asm_stats(str(fasta_file), buf, self.MONOMER, self.N)
        content = buf.getvalue()
        assert "CONTIG COUNT" in content
        assert "TOTAL LENGTH" in content

    def test_gz_fasta_input(self, gz_fasta_file, capsys):
        genome_asm_stats(str(gz_fasta_file), None, self.MONOMER, self.N)
        out = capsys.readouterr().out
        assert "CONTIG COUNT" in out
        assert "3" in out

    def test_telomere_counts(self, telomere_fasta_file, capsys):
        genome_asm_stats(str(telomere_fasta_file), None, self.MONOMER, self.N)
        out = capsys.readouterr().out
        # scaffold_t2t contributes to both fwd and t2t stats
        assert "T2T SCAFFOLDS" in out
        assert "TELOMERE FWD" in out
        assert "TELOMERE REV" in out


# ---------------------------------------------------------------------------
# run()
# ---------------------------------------------------------------------------

class TestAssessRun:
    def test_run_prints_to_stdout(self, assess_args, capsys):
        run(None, assess_args)
        out = capsys.readouterr().out
        assert "CONTIG COUNT" in out

    def test_run_writes_report_file(self, assess_args):
        run(None, assess_args)
        report_path = assess_args.report
        assert os.path.exists(report_path)
        content = open(report_path).read()
        assert "CONTIG COUNT" in content
        assert "TOTAL LENGTH" in content

    def test_run_no_report_arg(self, fasta_file, capsys):
        from argparse import Namespace
        args = Namespace(
            input=str(fasta_file),
            report=None,
            telomere_monomer="TAA[C]+",
            telomere_n_repeat=2,
            debug=False,
            pipe=True,
        )
        run(None, args)  # should not raise
        out = capsys.readouterr().out
        assert "CONTIG COUNT" in out
