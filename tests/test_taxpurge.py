"""Unit tests for aaftf/taxpurge.py (sourmash and its branchwater plugin are mocked)."""

import contextlib
import logging
import os
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.taxpurge import run

pytestmark = pytest.mark.unit

_TAX_HEADER = "query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank\n"
_FUNGUS = "Eukaryota;Ascomycota;Saccharomycetes"
_BACT = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria"


def _tax_row(contig, lineage, status="match", bp=5000):
    return f"{contig},{status},class,0.9,{lineage},md5,assembly.fasta,0.9,{bp},0.99\n"


@contextlib.contextmanager
def _chdir(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def _ids(path):
    return [line[1:].strip() for line in Path(path).read_text().splitlines() if line.startswith(">")]


def _fake_sourmash(rows, cmds, matched=True):
    """Return a run_cmd stand-in recording commands and writing what fastmultigather / tax genome would.

    ``rows`` are the classification rows, or a list of them per database (searched in order).
    """
    per_db = rows if rows and isinstance(rows[0], list) else None

    def _run_cmd(cmd, debug=False, cwd=None, **kw):
        cmds.append(cmd)
        if cmd[:3] == ["sourmash", "scripts", "fastmultigather"]:
            Path(cwd, cmd[cmd.index("-o") + 1]).write_text("intersect_bp,query_name\n" + ("3000,contig\n" if matched else ""))
        elif cmd[:3] == ["sourmash", "tax", "genome"]:
            prefix = cmd[cmd.index("-o") + 1]
            db_rows = per_db[int(prefix.split("_")[1])] if per_db else rows
            Path(cwd, f"{prefix}.classifications.csv").write_text(_TAX_HEADER + "".join(db_rows))

    return _run_cmd


def _run(tmp_path, rows, seqs=None, matched=True, **kw):
    """Run taxpurge on a small assembly with the given classification rows; return (kept ids, commands)."""
    seqs = seqs or {"fungal": "A" * 100, "bact": "C" * 100, "unknown": "G" * 100}
    asm = tmp_path / "asm.fasta"
    asm.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    (tmp_path / "db.sig.zip").write_text("")
    (tmp_path / "lineages.csv").write_text("ident,superkingdom\n")
    out = tmp_path / "out.fasta"
    cmds = []
    kw.setdefault("phylum", ["Ascomycota"])
    kw.setdefault("sourdb", [str(tmp_path / "db.sig.zip")])
    kw.setdefault("taxonomy_csv", [str(tmp_path / "lineages.csv")])
    with patch("aaftf.taxpurge.run_cmd", side_effect=_fake_sourmash(rows, cmds, matched)):
        # a relative work directory, as with `-w wd`: sourmash runs inside it
        with _chdir(tmp_path):
            run(input=str(asm), outfile=str(out), workdir="wd", **kw)
    return (_ids(out) if out.exists() else None), cmds


def test_drops_contigs_classified_outside_phylum(tmp_path):
    rows = [_tax_row("fungal", _FUNGUS), _tax_row("bact", _BACT)]
    kept, _ = _run(tmp_path, rows)
    assert kept == ["fungal", "unknown"]  # unclassified contigs are kept
    assert (tmp_path / "asm.sourmash-taxonomy.csv").read_text().startswith(_TAX_HEADER)


def test_below_threshold_is_kept(tmp_path):
    kept, _ = _run(tmp_path, [_tax_row("bact", _BACT, status="below_threshold")])
    assert kept == ["fungal", "bact", "unknown"]


def test_gtdb_rank_prefixes_are_optional_in_phylum(tmp_path):
    """GTDB lineages are p__-prefixed; -p works with or without the prefix."""
    rows = [_tax_row("fungal", _FUNGUS), _tax_row("bact", _BACT)]
    assert _run(tmp_path, rows, phylum=["Pseudomonadota"])[0] == ["bact", "unknown"]
    assert _run(tmp_path, rows, phylum=["p__Pseudomonadota"])[0] == ["bact", "unknown"]


def test_phylum_matches_names_exactly(tmp_path):
    rows = [_tax_row("fungal", _FUNGUS), _tax_row("bact", _BACT)]
    assert _run(tmp_path, rows, phylum=["Ascomy"])[0] == ["unknown"]


def test_header_with_description(tmp_path):
    """sourmash names each contig sketch after its full FASTA header; only the ID is matched."""
    kept, _ = _run(tmp_path, [_tax_row("bact some description", _BACT)])
    assert kept == ["fungal", "unknown"]


def test_sourmash_commands(tmp_path):
    _, (sketch, gather, tax) = _run(tmp_path, [], cpus=4)
    assert sketch == ["sourmash", "sketch", "dna", "-p", "k=31,scaled=1000", "--singleton", "-o", "contigs.sig.zip", "assembly.fasta"]
    assert gather == ["sourmash", "scripts", "fastmultigather", "contigs.sig.zip", str(tmp_path / "db.sig.zip"), "-k", "31", "-s", "1000", "-t", "2000", "-c", "4", "-o", "gather_0.csv"]
    assert tax == ["sourmash", "tax", "genome", "-g", "gather_0.csv", "-t", str(tmp_path / "lineages.csv"), "-F", "csv_summary", "-o", "tax_0"]


def test_databases_searched_one_at_a_time(tmp_path):
    """Each database gets its own fastmultigather run, so memory is that of one database."""
    (tmp_path / "db2.sig.zip").write_text("")
    _, cmds = _run(tmp_path, [], sourdb=[str(tmp_path / "db.sig.zip"), str(tmp_path / "db2.sig.zip")])
    gathers = [c for c in cmds if c[:3] == ["sourmash", "scripts", "fastmultigather"]]
    assert [(c[4], c[-1]) for c in gathers] == [(str(tmp_path / "db.sig.zip"), "gather_0.csv"), (str(tmp_path / "db2.sig.zip"), "gather_1.csv")]


def test_best_match_across_databases_wins(tmp_path):
    """A target contig partly matched by a (contaminated) bacterial genome keeps its fungal classification."""
    (tmp_path / "db2.sig.zip").write_text("")
    bacteria = [_tax_row("fungal", _BACT, bp=3000), _tax_row("bact", _BACT, bp=20000)]
    fungi = [_tax_row("fungal", _FUNGUS, bp=30000)]
    kept, _ = _run(tmp_path, [bacteria, fungi], sourdb=[str(tmp_path / "db.sig.zip"), str(tmp_path / "db2.sig.zip")])
    assert kept == ["fungal", "unknown"]
    saved = (tmp_path / "asm.sourmash-taxonomy.csv").read_text()
    assert "fungal,match" in saved and _FUNGUS in saved


def test_min_hashes_and_kmer_passed(tmp_path):
    _, (sketch, gather, _) = _run(tmp_path, [], min_hashes=5, kmer="21")
    assert "k=21,scaled=1000" in sketch
    assert gather[gather.index("-t") + 1] == "5000" and gather[gather.index("-k") + 1] == "21"


def _write_collection(path, scaled):
    """Write a minimal sourmash zip collection whose manifest gives ``scaled``."""
    import zipfile

    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("SOURMASH-MANIFEST.csv", f"# SOURMASH-MANIFEST-VERSION: 1.0\ninternal_location,md5,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\nsigs/a.sig.gz,x,31,DNA,0,{scaled},5,0,genome,f\n")


def test_database_searched_at_its_own_scaled(tmp_path):
    """A scaled=10000 collection is searched at -s 10000, with the bp threshold scaled with it."""
    _write_collection(tmp_path / "lite.sig.zip", 10000)
    _, (_, gather, _) = _run(tmp_path, [], sourdb=[str(tmp_path / "lite.sig.zip")])
    assert gather[gather.index("-s") + 1] == "10000" and gather[gather.index("-t") + 1] == "20000"


def test_no_match_keeps_everything_and_skips_tax(tmp_path):
    kept, cmds = _run(tmp_path, [], matched=False)
    assert kept == ["fungal", "bact", "unknown"]
    assert not any(c[:3] == ["sourmash", "tax", "genome"] for c in cmds)


def test_taxonomy_only_returns_early(tmp_path):
    """--just-show-taxonomy writes no filtered assembly, but saves the classification CSV and cleans up."""
    with patch("aaftf.taxpurge.cleanup_workdir") as cleanup:
        kept, _ = _run(tmp_path, [_tax_row("bact", _BACT)], taxonomy=True)
    assert kept is None
    assert (tmp_path / "asm.sourmash-taxonomy.csv").exists()
    cleanup.assert_called_once()


def _fake_require(names, hint=None):
    return [f"/db/{name}" for name in names]


def _fake_index(name):
    """find_db_index stand-in: the GTDB databases (published only as indexes) are installed."""
    return f"/db/{name}" if name.startswith("gtdb-") else None


def test_installed_database_names(tmp_path):
    """Each name gets its own scaled; lite and full databases sharing a lineage CSV pass it once."""
    with patch("aaftf.taxpurge.require_databases", side_effect=_fake_require), patch("aaftf.taxpurge.find_db_index", side_effect=_fake_index):
        _, cmds = _run(tmp_path, [], sourdb=["gbk-fungi", "gtdb-reps", "gtdb-reps-lite"], taxonomy_csv=None)
    gathers = [c for c in cmds if c[:3] == ["sourmash", "scripts", "fastmultigather"]]
    assert [(c[4], c[c.index("-s") + 1], c[c.index("-t") + 1]) for c in gathers] == [
        ("/db/gbk-fungi", "1000", "2000"),  # its sketch file (kept with --no-build): no index installed
        ("/db/gtdb-rs226-reps.k31.rocksdb", "1000", "2000"),
        ("/db/gtdb-rs226-reps.k31-sc10k.sig.rocksdb", "10000", "20000"),
    ]
    tax = next(c for c in cmds if c[:3] == ["sourmash", "tax", "genome"])
    assert tax[tax.index("-t") + 1 : tax.index("-F")] == ["/db/gbk-fungi-lineages", "/db/gtdb-reps-lineages"]


def test_human_hg38_lineage_written(tmp_path):
    """hg38-entire is published without a lineage CSV, so taxpurge writes one."""
    with patch("aaftf.taxpurge.require_databases", side_effect=_fake_require):
        _, cmds = _run(tmp_path, [], sourdb=["human-hg38"], taxonomy_csv=None)
    tax = next(c for c in cmds if c[:3] == ["sourmash", "tax", "genome"])
    lineage_csv = Path(tax[tax.index("-t") + 1])
    assert "hg38-entire,Eukaryota,Chordata" in lineage_csv.read_text()


def _install(db_dir, *names):
    """Put empty stand-ins for these DATABASES entries (sketch files, lineage CSVs) in ``db_dir``."""
    from aaftf.resources import DATABASES

    for name in names:
        (db_dir / DATABASES[name]["filename"]).write_text("")


def test_rocksdb_index_used_when_present(tmp_path, monkeypatch):
    """With its RocksDB index installed, a named database is searched from the index; its .zip is not needed."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    monkeypatch.setenv("AAFTF_DB", str(db_dir))
    _install(db_dir, "gtdb-reps-lineages")
    index = db_dir / "gtdb-rs226-reps.k31.rocksdb"
    index.mkdir()
    (index / "CURRENT").write_text("MANIFEST-000005\n")
    _, cmds = _run(tmp_path, [], sourdb=["gtdb-reps"], taxonomy_csv=None)
    gather = next(c for c in cmds if c[:3] == ["sourmash", "scripts", "fastmultigather"])
    assert gather[4] == str(db_dir.resolve() / "gtdb-rs226-reps.k31.rocksdb") and gather[gather.index("-s") + 1] == "1000"


def test_sketch_file_used_without_rocksdb_index(tmp_path, monkeypatch):
    """After `AAFTF database --no-build`, the sketch file is searched."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    monkeypatch.setenv("AAFTF_DB", str(db_dir))
    _install(db_dir, "gbk-fungi-lite", "gbk-euks-lineages")
    (db_dir / "genbank-euks-2025.01-fungi.k31.rocksdb").mkdir()  # incomplete: no CURRENT file
    _, cmds = _run(tmp_path, [], sourdb=["gbk-fungi-lite"], taxonomy_csv=None)
    gather = next(c for c in cmds if c[:3] == ["sourmash", "scripts", "fastmultigather"])
    assert gather[4] == str(db_dir.resolve() / "genbank-euks-2025.01-fungi.k31.sig.zip")


def test_published_index_database_missing_says_how_to_get_it(tmp_path, monkeypatch):
    monkeypatch.setenv("AAFTF_DB", str(tmp_path / "db"))
    with pytest.raises(FileNotFoundError, match="AAFTF database gtdb-reps-lite"):
        _run(tmp_path, [], sourdb=["gtdb-reps-lite"], taxonomy_csv=None)


def test_default_databases(tmp_path):
    with patch("aaftf.taxpurge.require_databases", side_effect=_fake_require) as req, patch("aaftf.taxpurge.find_db_index", side_effect=_fake_index):
        _run(tmp_path, [], sourdb=None, taxonomy_csv=None)
    assert [c.args[0] for c in req.call_args_list] == [["gbk-fungi-lite", "gbk-euks-lineages"], ["gtdb-reps-lineages"], ["human-hg38"]]


def test_database_path_without_lineages_raises(tmp_path):
    with pytest.raises(ValueError, match="--taxonomy_csv"):
        _run(tmp_path, [], taxonomy_csv=None)


def test_pipe_hides_hint(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    _run(tmp_path, [], pipe=True)
    assert "AAFTF covpurge" not in caplog.text
    _run(tmp_path, [], pipe=False)
    assert "AAFTF covpurge" in caplog.text
