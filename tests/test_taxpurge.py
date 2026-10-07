"""Unit tests for aaftf/taxpurge.py (sourmash is mocked)."""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.taxpurge import run

pytestmark = pytest.mark.unit

_HEADER = "ID,status,superkingdom,phylum,class,order,family,genus,species,strain\n"
_FUNGUS = "Eukaryota,Ascomycota,Sordariomycetes,,,,,\n"
_BACT = "Bacteria,Proteobacteria,Gammaproteobacteria,,,,,\n"
_NOMATCH = ",,,,,,,\n"


def _setup(tmp_path, seqs):
    asm = tmp_path / "asm.fasta"
    asm.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    db = tmp_path / "lca.json"
    db.write_text("{}")
    return asm, db


def _ids(path):
    return [line[1:].strip() for line in Path(path).read_text().splitlines() if line.startswith(">")]


def _run(tmp_path, classify, seqs=None, **kw):
    """Run taxpurge with the given classify CSV rows; return (kept ids, run_cmd mock, execute mock)."""
    seqs = seqs or {"fungal": "A" * 100, "bact": "C" * 100, "unknown": "G" * 100}
    asm, db = _setup(tmp_path, seqs)
    out = tmp_path / "out.fasta"

    def fake_execute(cmd, cwd=None, debug=False, quiet=False):
        yield from [_HEADER, *classify]

    with (
        patch("aaftf.taxpurge.run_cmd") as rc,
        patch("aaftf.taxpurge.execute", side_effect=fake_execute) as ex,
    ):
        kw.setdefault("phylum", ["Ascomycota"])
        run(input=str(asm), outfile=str(out), sourdb=str(db), workdir=str(tmp_path / "wd"), **kw)
    return (_ids(out) if out.exists() else None), rc, ex


def test_taxonomy_filter_keeps_phylum_and_nomatch(tmp_path):
    rows = ["fungal,found," + _FUNGUS, "bact,found," + _BACT, "unknown,nomatch," + _NOMATCH]
    kept, rc, ex = _run(tmp_path, rows)
    assert kept == ["fungal", "unknown"]
    rc.assert_called_once()
    assert rc.call_args.args[0][:2] == ["sourmash", "compute"]
    assert ex.call_args.args[0][:5] == ["sourmash", "lca", "classify", "--db", str(tmp_path / "lca.json")]
    assert (tmp_path / "asm.sourmash-taxonomy.csv").read_text().startswith(_HEADER)


def test_status_must_be_exactly_found(tmp_path):
    """Only status "found" counts as a classification; a status merely containing it is ignored."""
    kept, *_ = _run(tmp_path, ["fungal,found," + _FUNGUS, "bact,notfound," + _BACT])
    assert kept == ["fungal", "bact", "unknown"]


def test_phylum_matches_any_rank_exactly(tmp_path):
    # matching is exact membership in the lineage list, so a class name works and a prefix does not
    rows = ["fungal,found," + _FUNGUS, "bact,found," + _BACT]
    kept, *_ = _run(tmp_path, rows, phylum=["Gammaproteobacteria", "Ascomy"])
    assert kept == ["bact", "unknown"]


def test_taxonomy_only_returns_early(tmp_path):
    """--just-show-taxonomy writes no filtered assembly, but saves the classification CSV and cleans up."""
    with patch("aaftf.taxpurge.cleanup_workdir") as cleanup:
        kept, *_ = _run(tmp_path, ["bact,found," + _BACT], taxonomy=True)
    assert kept is None
    assert (tmp_path / "asm.sourmash-taxonomy.csv").read_text().startswith(_HEADER)
    cleanup.assert_called_once()  # the work directory (and its step log) is still cleaned up


def test_default_db_from_sourdb_type(tmp_path):
    asm, _ = _setup(tmp_path, {"a": "A" * 10})
    with (
        patch("aaftf.taxpurge.require_databases", return_value=[str(tmp_path / "gtdb.json")]) as req,
        patch("aaftf.taxpurge.run_cmd"),
        patch("aaftf.taxpurge.execute", return_value=iter([])) as ex,
    ):
        run(input=str(asm), outfile=str(tmp_path / "o.fa"), phylum=["X"], sourdb_type="gtdbrep", workdir=str(tmp_path / "wd"))
    assert req.call_args.args[0] == ["sm_gtdbrep"]
    assert str(tmp_path / "gtdb.json") in ex.call_args.args[0]


def test_pipe_hides_hint(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    _run(tmp_path, [], pipe=True)
    assert "AAFTF covpurge" not in caplog.text
    _run(tmp_path, [], pipe=False)
    assert "AAFTF covpurge" in caplog.text
