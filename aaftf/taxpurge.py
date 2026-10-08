"""Remove contigs classified outside the target taxa with sourmash.

Each contig is sketched, searched against sourmash databases with the branchwater plugin's
``fastmultigather``, and classified from those matches with ``sourmash tax genome``.
"""

import csv
import logging
import re
import shutil
import zipfile
from pathlib import Path
from typing import Any

from aaftf.resources import SOURMASH_DATABASES, TAXPURGE_DEFAULT_DATABASES
from aaftf.utility import cleanup_workdir, close_step_log, db_dirs, fasta_stats, filter_fasta, find_db_index, make_workdir, next_step_name, require_databases, run_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)

# contigs are sketched at this resolution (one hash kept per SCALED bp); each database is searched
# at its own, coarser or equal, resolution
SCALED = 1000

# lineage of the human-hg38 sketch, which is published without a lineage CSV
_HUMAN_LINEAGE = "ident,superkingdom,phylum,class,order,family,genus,species\nhg38-entire,Eukaryota,Chordata,Mammalia,Primates,Hominidae,Homo,Homo sapiens\n"

# GTDB lineages prefix each name with its rank ("p__Pseudomonadota"); dropped before matching --phylum
_RANK_PREFIX = re.compile(r"^[a-z]__")


def run(
    input: str,
    outfile: str,
    phylum: list[str],
    workdir: str | None = None,
    sourdb: list[str] | None = None,
    taxonomy_csv: list[str] | None = None,
    kmer: str = "31",
    min_hashes: int = 2,
    cpus: int = 1,
    taxonomy: bool = False,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Classify each contig with sourmash and drop those classified outside ``phylum``.

    Each contig is sketched (``sourmash sketch dna --singleton``, ``scaled=1000``), searched against
    each ``sourdb`` database in turn, at that database's resolution, with ``sourmash scripts
    fastmultigather`` (branchwater plugin), classified with ``sourmash tax genome``, and given the
    classification with the most matched sequence across databases. Classified contigs whose
    lineage contains none of ``phylum`` are dropped; unclassified contigs are kept. The
    classification CSV is copied next to ``input`` as ``<input stem>.sourmash-taxonomy.csv``.

    Args:
        input: Assembly FASTA to screen.
        outfile: Output FASTA of retained contigs.
        phylum: Taxon names (usually phyla) to keep; GTDB rank prefixes such as ``p__`` are optional.
        workdir: Working directory; a temporary one is created when None.
        sourdb: Databases to search: names in ``SOURMASH_DATABASES`` (installed by ``AAFTF
            database``), or paths of sourmash sketch collections or branchwater RocksDB indexes;
            default ``TAXPURGE_DEFAULT_DATABASES``.
        taxonomy_csv: Lineage CSV files for databases given as paths.
        kmer: k-mer size (as a string) of the sketches; must be one the databases contain.
        min_hashes: Minimum number of sketch hashes a contig must share with a database genome to
            match it, i.e. ``min_hashes`` x the database's ``scaled`` bp of estimated overlap;
            contigs with less stay unclassified.
        cpus: Threads for ``fastmultigather``.
        taxonomy: Only report the taxonomic classifications and return without filtering.
        debug: Print per-contig details and keep the work directory.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.

    Raises:
        ValueError: If no lineage CSV is available for the databases.
    """
    workdir, custom_workdir = make_workdir(workdir, "taxpurge")
    try:
        databases, lineages = _resolve_databases(sourdb or TAXPURGE_DEFAULT_DATABASES, taxonomy_csv or [], workdir)

        assembly_working = "assembly.fasta"
        shutil.copyfile(input, str(Path(workdir, assembly_working)))
        num_seqs, assembly_size = fasta_stats(str(Path(workdir, assembly_working)))
        logger.info(f"Assembly is {num_seqs:,} contigs and {assembly_size:,} bp")

        logger.info("Running sourmash to classify each contig")
        run_cmd(["sourmash", "sketch", "dna", "-p", f"k={kmer},scaled={SCALED}", "--singleton", "-o", "contigs.sig.zip", assembly_working], debug, cwd=workdir)
        # one database at a time, so memory is that of the largest database, not their sum; runs
        # in the work directory, since with sketch collections it also writes a CSV per contig
        gathers = []
        for index, (database, scaled) in enumerate(databases):
            gather = f"gather_{index}.csv"
            logger.info(f"Searching {database} (scaled={scaled}: contigs need {min_hashes * scaled:,} bp matching a genome)")
            run_cmd(["sourmash", "scripts", "fastmultigather", "contigs.sig.zip", database, "-k", kmer, "-s", str(scaled), "-t", str(min_hashes * scaled), "-c", str(cpus), "-o", gather], debug, cwd=workdir)
            gathers.append(gather)
        contig_taxonomy, classification_csv = _classify(workdir, gathers, lineages, debug)

        phyla = {}
        for lineage in contig_taxonomy.values():
            name = ";".join(lineage[:2])
            phyla[name] = phyla.get(name, 0) + 1
        summary = "\n".join(f"{count:>6,}  {name}" for name, count in sorted(phyla.items(), key=lambda item: -item[1]))
        logger.info(f"Classified {len(contig_taxonomy):,} of {num_seqs:,} contigs:\n{summary}")
        if classification_csv:
            _save_taxonomy_csv(classification_csv, input)
        if taxonomy:  # --just-show-taxonomy only reports the classifications
            cleanup_workdir(workdir, debug, custom_workdir)
            return

        keep = {_RANK_PREFIX.sub("", name) for name in phylum}
        tax_to_drop = {contig for contig, lineage in contig_taxonomy.items() if not keep.intersection(lineage)}
        if debug:
            for contig in sorted(tax_to_drop):
                logger.info(f"{contig}\t{';'.join(contig_taxonomy[contig])}")

        logger.info(f"Dropping {len(tax_to_drop):,} contigs classified outside {', '.join(phylum)}")
        num_seqs, assembly_size = filter_fasta(str(Path(workdir, assembly_working)), outfile, lambda seq_id: seq_id not in tax_to_drop)
        logger.info(f"Taxonomy-purged assembly is {num_seqs:,} contigs and {assembly_size:,} bp")
        cleanup_workdir(workdir, debug, custom_workdir)

        if not pipe:
            next_out = next_step_name(outfile, ".covpurge.fasta")
            logger.info(f"Your next command might be:\nAAFTF covpurge -i {outfile} -o {next_out} -1 READ1 -2 READ2")
    finally:
        close_step_log()


def _resolve_databases(sourdb: list[str], taxonomy_csv: list[str], workdir: str) -> tuple[list[tuple[str, int]], list[str]]:
    """Return the ((sketch collection, scaled) pairs, lineage CSVs) for the ``sourdb`` names and paths.

    Args:
        sourdb: Names in ``SOURMASH_DATABASES`` or paths of sketch collections / RocksDB indexes.
        taxonomy_csv: Lineage CSVs for the databases given as paths.
        workdir: Work directory, where the human-hg38 lineage CSV is written.

    Raises:
        ValueError: If no lineage CSV is available.
    """
    databases: list[tuple[str, int]] = []
    lineages = [str(Path(f).resolve()) for f in taxonomy_csv]
    for item in sourdb:
        if item in SOURMASH_DATABASES:
            entry = SOURMASH_DATABASES[item]
            # its RocksDB index (what `AAFTF database` installs) is searched from disk; the sketch
            # file, kept by `AAFTF database --no-build`, is loaded into memory
            index = find_db_index(entry["rocksdb"])
            if not index and "filename" not in entry:  # only published as an index
                raise FileNotFoundError(f"missing database: {item} ({entry['rocksdb']}; searched {', '.join(str(d) for d in db_dirs())})\nDownload it first: AAFTF database {item}")
            names = ([] if index else [item]) + ([entry["lineages"]] if entry["lineages"] else [])
            paths = require_databases(names, hint=f"(AAFTF database {item} downloads both) or pass --sourdb PATH with --taxonomy_csv")
            if index:
                logger.info(f"Using the RocksDB index {index} for {item}")
                databases.append((index, entry["scaled"]))
            else:
                databases.append((paths.pop(0), entry["scaled"]))
            if entry["lineages"]:
                lineages.append(paths[0])
            else:  # human-hg38 is published without one
                human = Path(workdir, "human-hg38.lineages.csv")
                human.write_text(_HUMAN_LINEAGE)
                lineages.append(str(human.resolve()))
        else:
            path = str(Path(item).resolve())
            databases.append((path, _collection_scaled(path)))
    lineages = list(dict.fromkeys(lineages))  # databases may share a lineage CSV
    if not lineages:
        raise ValueError("no lineage CSV for the sourmash databases; pass --taxonomy_csv with --sourdb PATH")
    return databases, lineages


def _collection_scaled(path: str) -> int:
    """Return the coarsest ``scaled`` of a sketch collection (from its manifest), or ``SCALED``.

    A RocksDB index or a collection without a manifest is assumed to be at ``SCALED``.
    """
    try:
        with zipfile.ZipFile(path) as zf:
            manifest = zf.read("SOURMASH-MANIFEST.csv").decode().splitlines()
        return max(int(row["scaled"]) for row in csv.DictReader(line for line in manifest if not line.startswith("#")))
    except (OSError, KeyError, ValueError, zipfile.BadZipFile):
        return SCALED


def _classify(workdir: str, gathers: list[str], lineages: list[str], debug: bool) -> tuple[dict[str, list[str]], str | None]:
    """Classify the contigs from each database's ``fastmultigather`` results with ``sourmash tax genome``.

    Each database is classified separately; a contig matched in several databases takes the
    classification with the most matched sequence (``bp_match_at_rank``), so a contig of the
    target organism is not called a contaminant because a contaminant genome contains a little of
    it.

    Args:
        workdir: Work directory holding the ``gathers`` CSVs.
        gathers: ``fastmultigather`` result CSVs, one per database.
        lineages: Lineage CSVs for all databases.
        debug: Show sourmash output.

    Returns:
        ``{contig id: lineage names, rank prefixes removed}`` for the contigs classified (status
        ``match``), and the CSV of the chosen classifications (None when no contig matched).
    """
    best: dict[str, dict[str, str]] = {}
    for index, gather in enumerate(gathers):
        gather_path = Path(workdir, gather)
        if not gather_path.is_file() or sum(1 for _ in gather_path.open()) < 2:
            continue
        run_cmd(["sourmash", "tax", "genome", "-g", gather, "-t", *lineages, "-F", "csv_summary", "-o", f"tax_{index}"], debug, cwd=workdir)
        with open(Path(workdir, f"tax_{index}.classifications.csv"), newline="") as fh:
            for row in csv.DictReader(fh):
                if row["status"] != "match" or not row["lineage"]:
                    continue
                contig = row["query_name"].split()[0]
                if contig not in best or float(row["bp_match_at_rank"]) > float(best[contig]["bp_match_at_rank"]):
                    best[contig] = row
    if not best:
        logger.info("No contig matched the sourmash databases")
        return {}, None
    classification_csv = str(Path(workdir, "classifications.csv"))
    with open(classification_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(next(iter(best.values()))))
        writer.writeheader()
        writer.writerows(best.values())
    contig_taxonomy = {contig: [_RANK_PREFIX.sub("", name) for name in row["lineage"].split(";") if name] for contig, row in best.items()}
    return contig_taxonomy, classification_csv


def _save_taxonomy_csv(classification_csv: str, input: str) -> None:
    """Copy the classification CSV next to the input as ``<input stem>.sourmash-taxonomy.csv``.

    Args:
        classification_csv: The ``sourmash tax genome`` classification CSV in the work directory.
        input: The input assembly; the copy goes in its directory, named after it.
    """
    baseinput = Path(input).name
    if "." in baseinput:
        baseinput = baseinput.rsplit(".", 1)[0]
    shutil.copy(classification_csv, str(Path(Path(input).parent, baseinput + ".sourmash-taxonomy.csv")))
