"""Unit tests for AAFTF/utility.py.

All tests are pure Python — no external bioinformatics tools required.
"""

import gzip

import pytest

from AAFTF.utility import (
    RevComp,
    SafeRemove,
    calcN50,
    checkfile,
    countfasta,
    countfastq,
    fastastats,
    line_count,
    myround,
    softwrap,
    which,
)

from tests.conftest import SMALL_FASTA, make_fastq_text

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# myround
# ---------------------------------------------------------------------------

class TestMyround:
    def test_round_up_to_base10(self):
        assert myround(15, 10) == 20

    def test_round_down_to_base10(self):
        assert myround(14, 10) == 10

    def test_exact_multiple_unchanged(self):
        assert myround(20, 10) == 20

    def test_base5(self):
        assert myround(12, 5) == 10

    def test_base5_round_up(self):
        assert myround(13, 5) == 15

    def test_zero(self):
        assert myround(0, 10) == 0


# ---------------------------------------------------------------------------
# calcN50
# ---------------------------------------------------------------------------

class TestCalcN50:
    """calcN50(lengths) uses the default num=0.5 for N50."""

    # lengths [20, 80, 150] → total=250, 50 %=125
    # cumsum from largest: 150 ≥ 125 → N50=150
    def test_n50_basic(self):
        assert calcN50([20, 80, 150]) == 150

    # 90 %=225 → cumsum 150 < 225, 150+80=230 ≥ 225 → N90=80
    def test_n90(self):
        assert calcN50([20, 80, 150], num=0.9) == 80

    def test_single_contig(self):
        assert calcN50([500]) == 500

    def test_equal_contigs(self):
        # [100]*4 → total=400, 50 %=200; cumsum: 100 < 200, 200 ≥ 200 → N50=100
        assert calcN50([100, 100, 100, 100]) == 100

    def test_two_contigs(self):
        # [300, 100] → total=400, 50%=200; cumsum: 300 ≥ 200 → N50=300
        assert calcN50([300, 100]) == 300

    def test_input_order_irrelevant(self):
        # Function sorts internally
        assert calcN50([150, 20, 80]) == calcN50([20, 80, 150])


# ---------------------------------------------------------------------------
# RevComp
# ---------------------------------------------------------------------------

class TestRevComp:
    def test_simple(self):
        assert RevComp("ATCG") == "CGAT"

    def test_complement_only(self):
        assert RevComp("AAAA") == "TTTT"
        assert RevComp("CCCC") == "GGGG"

    def test_palindrome(self):
        assert RevComp("AATTAATT") == "AATTAATT"

    def test_uppercase_conversion(self):
        # RevComp uppercases input before processing
        assert RevComp("atcg") == "CGAT"

    def test_longer_sequence(self):
        assert RevComp("ATCGATCG") == "CGATCGAT"

    def test_all_bases(self):
        # A↔T, C↔G
        assert RevComp("ACGT") == "ACGT"


# ---------------------------------------------------------------------------
# softwrap
# ---------------------------------------------------------------------------

class TestSoftwrap:
    def test_short_string_unchanged(self):
        assert softwrap("ATCG", 80) == "ATCG"

    def test_wraps_at_specified_width(self):
        result = softwrap("A" * 20, 10)
        lines = result.split("\n")
        assert len(lines) == 2
        assert all(len(line) == 10 for line in lines)

    def test_exact_width_no_wrap(self):
        assert softwrap("A" * 10, 10) == "A" * 10

    def test_one_over_width_wraps(self):
        result = softwrap("A" * 11, 10)
        assert "\n" in result
        assert result.split("\n")[0] == "A" * 10

    def test_empty_string(self):
        assert softwrap("", 80) == ""


# ---------------------------------------------------------------------------
# checkfile
# ---------------------------------------------------------------------------

class TestCheckfile:
    def test_existing_nonempty_file_is_true(self, fasta_file):
        assert checkfile(str(fasta_file)) is True

    def test_empty_file_is_false(self, empty_file):
        assert checkfile(str(empty_file)) is False

    def test_missing_file_is_false(self, tmp_path):
        assert checkfile(str(tmp_path / "nonexistent.fa")) is False


# ---------------------------------------------------------------------------
# line_count
# ---------------------------------------------------------------------------

class TestLineCount:
    def test_known_line_count(self, line_file):
        assert line_count(str(line_file)) == 5

    def test_fasta_line_count(self, fasta_file):
        # SMALL_FASTA has 3 header lines + 3 sequence lines = 6
        assert line_count(str(fasta_file)) == 6


# ---------------------------------------------------------------------------
# countfasta / fastastats
# ---------------------------------------------------------------------------

class TestCountFasta:
    def test_count_three_sequences(self, fasta_file):
        assert countfasta(str(fasta_file)) == 3

    def test_fastastats_count(self, fasta_file):
        count, total_len = fastastats(str(fasta_file))
        assert count == 3

    def test_fastastats_total_length(self, fasta_file):
        # SEQ1=20, SEQ2=80, SEQ3=150 → 250
        count, total_len = fastastats(str(fasta_file))
        assert total_len == 250


# ---------------------------------------------------------------------------
# countfastq
# ---------------------------------------------------------------------------

class TestCountFastq:
    def test_plain_fastq(self, fastq_file):
        assert countfastq(str(fastq_file)) == 10

    def test_gzip_fastq(self, gz_fastq_file):
        assert countfastq(str(gz_fastq_file)) == 10

    def test_single_read(self, tmp_path):
        p = tmp_path / "one.fastq"
        p.write_text(make_fastq_text(1))
        assert countfastq(str(p)) == 1


# ---------------------------------------------------------------------------
# which
# ---------------------------------------------------------------------------

class TestWhich:
    def test_python_found(self):
        # python3 or python should be in PATH in any test environment
        result = which("python3") or which("python")
        assert result is not None

    def test_nonexistent_tool_returns_none(self):
        assert which("definitely_not_a_real_tool_xyzzy_12345") is None


# ---------------------------------------------------------------------------
# SafeRemove
# ---------------------------------------------------------------------------

class TestSafeRemove:
    def test_remove_file(self, tmp_path):
        p = tmp_path / "todelete.txt"
        p.write_text("bye")
        SafeRemove(str(p))
        assert not p.exists()

    def test_remove_directory(self, tmp_path):
        d = tmp_path / "subdir"
        d.mkdir()
        (d / "file.txt").write_text("hi")
        SafeRemove(str(d))
        assert not d.exists()

    def test_remove_nonexistent_is_noop(self, tmp_path):
        # Should not raise
        SafeRemove(str(tmp_path / "ghost"))
