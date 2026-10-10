"""Unit tests for aaftf/vecscreen.py (blastn/makeblastdb are mocked)."""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf import vecscreen

pytestmark = pytest.mark.unit


def _read_fasta(path):
    """Return {id: seq} for a FASTA file."""
    recs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = ""
        elif name:
            recs[name] += line.strip()
    return recs


def _row(qid, pident, length, qstart, qend):
    return [qid, "hit_" + qid, str(pident), str(length), "0", "0", str(qstart), str(qend), "1", "50", "1e-10", "100"]


# --- _screen_euk_prok_contamination -------------------------------------------------------------


def test_euk_prok_screen_filters_and_sorts():
    rows = {
        "CONTAM_EUKS": [
            _row("c1", 98.0, 50, 300, 100),  # qualifies (98/50), reversed coords -> sorted
            _row("c2", 97.9, 99, 1, 99),  # fails: <98 at <100bp
            _row("c3", 94.0, 100, 10, 110),  # qualifies (94/100)
        ],
        "CONTAM_PROKS": [
            _row("c1", 90.0, 200, 500, 700),  # qualifies (90/200), same contig grouped
            _row("c4", 93.9, 199, 1, 199),  # fails
        ],
    }

    def fake(infile, workdir, prefix, contam, cpus, pid):
        return rows[contam]

    with patch("aaftf.vecscreen._run_blastn_screen", side_effect=fake):
        regions = vecscreen._screen_euk_prok_contamination("in.fa", "wd", "p", 2, "90.0")
    assert regions == {
        "c1": [(100, 300, "CONTAM_EUKS", "hit_c1", 98.0), (500, 700, "CONTAM_PROKS", "hit_c1", 90.0)],
        "c3": [(10, 110, "CONTAM_EUKS", "hit_c3", 94.0)],
    }


def test_write_euk_cleaned_splits_regions(tmp_path):
    infile = tmp_path / "in.fa"
    seq = "A" * 10 + "C" * 5 + "G" * 10
    infile.write_text(f">c1\n{seq}\n>c2\nTTTT\n")
    assert vecscreen._write_euk_cleaned(str(infile), {}, str(tmp_path), "p") == str(infile)
    out = vecscreen._write_euk_cleaned(str(infile), {"c1": [(11, 15, "E", "h", 99.0)]}, str(tmp_path), "p")
    assert _read_fasta(out) == {"split0_c1": "A" * 10, "split1_c1": "G" * 10, "c2": "TTTT"}


def test_screen_mitochondria_flags_long_hits():
    to_remove = {}
    rows = [_row("m1", 99, 120, 1, 120), _row("m2", 99, 119, 1, 119)]
    with patch("aaftf.vecscreen._run_blastn_screen", return_value=rows) as blast:
        hits = vecscreen._screen_mitochondria("q.fa", "wd", "p", 1, to_remove)
    assert hits == ["m1"]
    assert to_remove == {"m1": ("MitoScreen", "hit_m1", 99.0)}
    assert blast.call_args.args[3:] == ("MITO", 1, vecscreen.BLAST_PERCENT_ID_MITO_MATCH)


# --- _classify_vector_hits ----------------------------------------------------------------------


def _vec_row(qid, qstart, qend, score, qlen=1000):
    return "\t".join([qid, "vec", "100", "30", "0", "0", str(qstart), str(qend), "1", "30", "1", "50", str(score), str(qlen)])


@pytest.mark.parametrize(
    ("qstart", "qend", "score", "high", "low"),
    [
        (1, 30, 24, True, True),  # terminal strong
        (1, 30, 19, True, False),  # terminal moderate
        (1, 30, 18, False, False),  # terminal weak
        (975, 1000, 24, True, True),  # 3' terminal strong
        (500, 530, 30, True, True),  # internal strong
        (500, 530, 25, True, False),  # internal moderate
        (500, 530, 24, False, False),  # internal weak (would be strong if terminal)
        (26, 60, 24, False, False),  # 26 > 25 -> internal, so 24 is weak
    ],
)
def test_classify_thresholds(tmp_path, qstart, qend, score, high, low):
    tab = tmp_path / "r.tab"
    tab.write_text(_vec_row("c1", qstart, qend, score) + "\n")
    for stringency, expected in (("high", high), ("low", low)):
        hits, n = vecscreen._classify_vector_hits(str(tab), stringency, {})
        assert n == int(expected)
        assert ("c1" in hits) == expected


def test_classify_positions_and_skip_removed(tmp_path):
    tab = tmp_path / "r.tab"
    tab.write_text("\n".join([_vec_row("c1", 40, 5, 30), _vec_row("c1", 990, 960, 30), _vec_row("c1", 500, 540, 30), _vec_row("gone", 1, 30, 40)]) + "\n")
    hits, n = vecscreen._classify_vector_hits(str(tab), "high", {"gone": ("MitoScreen", "x", 99.0)})
    assert n == 3
    assert hits == {
        "c1": [
            ("vec", 1000, [5, 40], 30, True, "5"),
            ("vec", 1000, [960, 990], 30, True, "3"),
            ("vec", 1000, [500, 540], 30, False, None),
        ]
    }


# --- _write_trimmed_and_split / _group ----------------------------------------------------------


def test_group_drops_short_tail():
    assert list(vecscreen._group([1, 2, 3, 4, 5], 2)) == [(1, 2), (3, 4)]


def test_terminal_trim_and_internal_split(tmp_path):
    fa = tmp_path / "in.fa"
    seq_t = "N" * 20 + "A" * 300 + "T" * 30  # 350bp: 5' vector 1..20, 3' vector 321..350
    seq_s = "A" * 250 + "G" * 50 + "C" * 250  # internal vector 251..300
    fa.write_text(f">term\n{seq_t}\n>split\n{seq_s}\n>short\n{'A' * 199}\n>ok\n{'A' * 200}\n")
    hits = {
        "term": [("v", 350, [1, 20], 30, True, "5"), ("v", 350, [321, 350], 30, True, "3")],
        "split": [("v", 550, [251, 300], 30, False, None), ("v", 550, [251, 300], 30, False, None)],
    }
    out = tmp_path / "out.fsa"
    vecscreen._write_trimmed_and_split(str(fa), hits, str(out))
    assert _read_fasta(out) == {"term": "A" * 300, "split1_split": "A" * 250, "split2_split": "C" * 250, "ok": "A" * 200}


def test_split_pieces_shorter_than_200_dropped(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">c\n" + "A" * 100 + "G" * 50 + "C" * 250 + "\n")
    out = tmp_path / "out.fsa"
    vecscreen._write_trimmed_and_split(str(fa), {"c": [("v", 400, [101, 150], 30, False, None)]}, str(out))
    assert _read_fasta(out) == {"split2_c": "C" * 250}


def test_parse_clean_blastn(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">c\n" + "A" * 300 + "\n")
    tab = tmp_path / "r.tab"
    tab.write_text(_vec_row("c", 1, 30, 30, qlen=300) + "\n")
    n, cleaned = vecscreen._parse_clean_blastn(str(fa), str(tmp_path / "pre"), str(tab), "high")
    assert (n, cleaned) == (1, str(tmp_path / "pre") + ".clean.fsa")
    assert _read_fasta(cleaned) == {"c": "A" * 270}


# --- _write_final_outputs -----------------------------------------------------------------------


def test_write_final_outputs(tmp_path):
    vec = tmp_path / "vec.fa"
    vec.write_text(">keep\nAAAA\n>mito\nCCCC\n>contam\nGGGG\n")
    out, mito = tmp_path / "clean.fa", tmp_path / "mito.fa"
    vecscreen._write_final_outputs(str(vec), {"mito": ("MitoScreen", "h", 99.0), "contam": ("X", "h", 99.0)}, ["mito"], str(out), str(mito))
    assert _read_fasta(out) == {"keep": "AAAA"}
    assert _read_fasta(mito) == {"mito": "CCCC"}


# --- _make_blastdb / _run_vecscreen_rounds ------------------------------------------------------


def test_make_blastdb_skips_up_to_date_index(tmp_path):
    fa = tmp_path / "db.fa"
    fa.write_text(">a\nA\n")
    with patch("aaftf.vecscreen.call") as call:
        vecscreen._make_blastdb("nucl", str(fa), str(tmp_path / "db"))
        assert call.call_args.args[0][:3] == ["makeblastdb", "-dbtype", "nucl"]
        Path(tmp_path / "db.nin").write_text("")
        call.reset_mock()
        vecscreen._make_blastdb("nucl", str(fa), str(tmp_path / "db"))
        call.assert_not_called()
        vecscreen._make_blastdb("prot", str(fa), str(tmp_path / "db"))
        call.assert_called_once()


def test_vecscreen_rounds_repeat_until_clean(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">c\n" + "A" * 300 + "\n")
    reports = iter([_vec_row("c", 1, 30, 30, qlen=300) + "\n", ""])
    queries = []

    def fake_call(cmd, **kw):
        queries.append(cmd[cmd.index("-query") + 1])
        Path(cmd[cmd.index("-out") + 1]).write_text(next(reports))

    with patch("aaftf.vecscreen.call", side_effect=fake_call):
        final = vecscreen._run_vecscreen_rounds(str(fa), str(tmp_path), "p", 1, "high", {})
    assert queries == [str(fa), str(tmp_path / "p.r0.clean.fsa")]
    assert final == str(tmp_path / "p.r1.clean.fsa")
    assert _read_fasta(final) == {"c": "A" * 270}


# --- run() --------------------------------------------------------------------------------------


@pytest.mark.parametrize("pipe", [False, True])
def test_run_end_to_end(tmp_path, monkeypatch, caplog, pipe):
    monkeypatch.chdir(tmp_path)
    infile = tmp_path / "asm.fasta"
    infile.write_text(">keep\n" + "A" * 300 + "\n>mito\n" + "C" * 300 + "\n")
    srcs = []
    for name in ("univec", "euks", "proks", "mitodb"):
        p = tmp_path / f"{name}.fa"
        p.write_text(">x\nACGT\n")
        srcs.append(str(p))
    wd = tmp_path / "wd"
    blast_calls = []

    def fake_call(cmd, **kw):
        if cmd[0] == "makeblastdb":
            return 0
        out = cmd[cmd.index("-out") + 1]
        blast_calls.append(Path(out).name)
        rows = "mito\tmt\t99.5\t300\t0\t0\t1\t300\t1\t300\t0\t500\n" if "MITO" in out else ""
        Path(out).write_text(rows)
        return 0

    caplog.set_level(logging.INFO)
    (tmp_path / "sub").mkdir()
    with patch("aaftf.vecscreen.require_databases", return_value=srcs), patch("aaftf.vecscreen.call", side_effect=fake_call):
        vecscreen.run(infile=str(infile), outfile="sub/asm.vecscreen.fasta", workdir=str(wd), pipe=pipe)
    assert blast_calls == ["CONTAM_EUKS.asm.vecscreen.blastn", "CONTAM_PROKS.asm.vecscreen.blastn", "MITO.asm.vecscreen.blastn", "asm.vecscreen.r0.vecscreen.tab"]
    # both outputs go next to -o/--outfile, not into the current directory
    assert _read_fasta(tmp_path / "sub" / "asm.vecscreen.fasta") == {"keep": "A" * 300}
    assert _read_fasta(tmp_path / "sub" / "asm.vecscreen.mitochondria.fasta") == {"mito": "C" * 300}
    assert not (tmp_path / "asm.vecscreen.fasta").exists()
    assert ("AAFTF taxpurge" in caplog.text) is (not pipe)
