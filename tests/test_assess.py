"""Unit tests for aaftf/assess.py.

Tests cover the pure-Python helper functions (make_regex_revcomp, find_telomere)
and the full genome_asm_stats / run() pipeline against known small
FASTA files.  No external tools are required.
"""

import inspect
import io
from pathlib import Path
from unittest.mock import patch

import pytest
from Bio.Seq import Seq

from aaftf.assess import find_telomere, genome_asm_stats, make_regex_revcomp, run
from tests.conftest import (
    TELO_FWD_ONLY,
    TELO_NONE,
    TELO_REV_ONLY,
    TELO_T2T,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# make_regex_revcomp  (regex-aware)
# ---------------------------------------------------------------------------


class TestRevcomp:
    """make_regex_revcomp handles regex bracket metacharacters in the monomer."""

    def test_simple_dna(self):
        assert make_regex_revcomp("ATCG") == "CGAT"

    def test_all_complement(self):
        assert make_regex_revcomp("AAAA") == "TTTT"
        assert make_regex_revcomp("CCCC") == "GGGG"

    def test_regex_monomer(self):
        # revcomp of "TAA[C]+" should be "[G]+TTA"
        rc = make_regex_revcomp("TAA[C]+")
        assert rc == "[G]+TTA"

    def test_symmetric(self):
        # make_regex_revcomp(make_regex_revcomp(seq)) should equal seq for plain DNA
        seq = "ATCGATCG"
        assert make_regex_revcomp(make_regex_revcomp(seq)) == seq


# ---------------------------------------------------------------------------
# find_telomere
# ---------------------------------------------------------------------------


class TestFindTelomere:
    MONOMER = "TAA[C]+"
    N = 2

    def _check(self, seq_str):
        from Bio.Seq import Seq

        return find_telomere(Seq(seq_str), self.MONOMER, self.N)

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
        run(**vars(assess_args))
        out = capsys.readouterr().out
        assert "CONTIG COUNT" in out

    def test_run_writes_report_file(self, assess_args):
        run(**vars(assess_args))
        report_path = assess_args.report
        assert Path(report_path).exists()
        content = Path(report_path).read_text()
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
        )
        run(**vars(args))  # should not raise
        out = capsys.readouterr().out
        assert "CONTIG COUNT" in out


class TestAssessRunErrors:
    def test_missing_input_raises_before_writing_report(self, tmp_path):
        report = tmp_path / "report.txt"
        with pytest.raises(FileNotFoundError, match="assembly file not found"):
            run(input=str(tmp_path / "missing.fasta"), report=str(report))
        assert not report.exists()

    @pytest.mark.parametrize("content", ["", ">empty\n"])
    def test_empty_assembly_raises(self, tmp_path, content):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text(content)
        with pytest.raises(ValueError, match="contains no sequence"):
            run(input=str(fasta))


# ---------------------------------------------------------------------------
# find_telomere with the real CLI default monomer
# ---------------------------------------------------------------------------


class TestFindTelomereDefaultMonomer:
    """Telomere calls with run()'s default monomer and repeat count."""

    MONOMER = inspect.signature(run).parameters["telomere_monomer"].default
    N = inspect.signature(run).parameters["telomere_n_repeat"].default
    FILLER = "GATTGCATGCAGTCAGTGCA" * 30  # 600 bp, no TAAC

    def test_default_repeat_count(self):
        assert self.N == 2

    def test_forward_telomere(self):
        seq = "TAACCC" * 4 + self.FILLER
        assert find_telomere(Seq(seq), self.MONOMER, self.N) == (True, False)

    def test_reverse_telomere(self):
        seq = str(Seq("TAACCC" * 4 + self.FILLER).reverse_complement())
        assert find_telomere(Seq(seq), self.MONOMER, self.N) == (False, True)

    def test_single_copy_not_called(self):
        seq = "TAACCC" + self.FILLER
        assert find_telomere(Seq(seq), self.MONOMER, self.N) == (False, False)

    def test_taac_with_one_c_not_called(self):
        seq = "TAAC" * 6 + self.FILLER
        rc = str(Seq(seq).reverse_complement())
        assert find_telomere(Seq(seq), self.MONOMER, self.N) == (False, False)
        assert find_telomere(Seq(rc), self.MONOMER, self.N) == (False, False)


class TestAssessRunReportHandle:
    def test_report_closed_when_stats_raise(self, fasta_file, tmp_path):
        report = tmp_path / "report.txt"
        handles = []

        def boom(fasta, handle, *args):
            handles.append(handle)
            raise RuntimeError("fail")

        with patch("aaftf.assess.genome_asm_stats", side_effect=boom):
            with pytest.raises(RuntimeError):
                run(input=str(fasta_file), report=str(report))
        assert handles[0].closed
        assert report.read_text() == ""


class TestSoftMaskReport:
    """BASES MASKED / PERCENT MASKED are reported only when the assembly has soft-masked bases."""

    def _report(self, tmp_path, seq):
        fasta = tmp_path / "a.fasta"
        fasta.write_text(f">c\n{seq}\n")
        report = tmp_path / "r.txt"
        run(input=str(fasta), report=str(report))
        return report.read_text()

    def test_no_lowercase_no_mask_lines(self, tmp_path):
        assert "MASKED" not in self._report(tmp_path, "ACGT" * 50)

    def test_some_lowercase_reported(self, tmp_path):
        text = self._report(tmp_path, "ACGT" * 45 + "acgt" * 5)
        assert "BASES MASKED  =  20" in text and "PERCENT MASKED  =  10.00" in text
