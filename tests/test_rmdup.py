"""Unit tests for aaftf/rmdup.py (minimap2 is stubbed; no external tools needed)."""

from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.rmdup import _is_duplicate, _write_query_and_reference, run
from aaftf.utility import PafHit

pytestmark = pytest.mark.unit


def _hit(query_len, matches, aln_len):
    return PafHit("q", query_len, 0, aln_len, "+", "t", 10_000, 0, aln_len, matches, aln_len, 60)


class TestWriteQueryAndReference:
    def test_exact_id_match_not_substring(self, tmp_path):
        # ctg1 is a substring of ctg10; only ctg10 may go to the query file
        seqs = {"ctg1": "AAAA", "ctg10": "CCCC", "ctg2": "GGGG"}
        query, reference = _write_query_and_reference(seqs, "ctg10", ["ctg1", "ctg2"], str(tmp_path), "p")
        assert Path(query).read_text() == ">ctg10\nCCCC\n"
        assert Path(reference).read_text() == ">ctg1\nAAAA\n>ctg2\nGGGG\n"


class TestIsDuplicate:
    def _check(self, hits, percent_id=95, percent_cov=95):
        with patch("aaftf.rmdup.paf_hits", return_value=iter(hits)):
            return _is_duplicate("q.fa", "r.fa", "q", 1, percent_id, percent_cov, False)

    def test_high_identity_and_coverage_is_duplicate(self):
        assert self._check([_hit(query_len=1000, matches=990, aln_len=1000)])

    def test_low_coverage_is_not(self):
        assert not self._check([_hit(query_len=1000, matches=500, aln_len=500)])

    def test_low_identity_is_not(self):
        assert not self._check([_hit(query_len=1000, matches=800, aln_len=1000)])

    def test_no_hits_is_not(self):
        assert not self._check([])


class TestRun:
    def test_drops_short_and_duplicate_contigs(self, tmp_path):
        fasta = tmp_path / "asm.fa"
        # ctg1 is short (< minlen); ctg10 will be reported duplicated; the rest are longest
        fasta.write_text(">ctg1\n" + "A" * 50 + "\n>ctg10 desc\n" + "C" * 600 + "\n>big1\n" + "G" * 5000 + "\n>big2\n" + "T" * 6000 + "\n")
        out = tmp_path / "out.fa"
        checked = []

        def _fake_is_duplicate(query, reference, name, *args):
            checked.append(name)
            return name == "ctg10"

        with patch("aaftf.rmdup._is_duplicate", side_effect=_fake_is_duplicate):
            run(input=str(fasta), out=str(out), workdir=str(tmp_path / "wd"), minlen=100)

        headers = [line for line in out.read_text().splitlines() if line.startswith(">")]
        assert headers == [">big1", ">big2"]
        # N75 is 5000 (big1's length): contigs up to and including N75 are checked, longer ones are not
        assert checked == ["ctg10", "big1"]
