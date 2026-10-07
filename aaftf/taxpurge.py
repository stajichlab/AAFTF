"""Remove contigs classified outside the target taxa with sourmash."""

import logging
import shutil
from pathlib import Path
from typing import Any

from aaftf.utility import check_file, cleanup_workdir, close_step_log, execute, fasta_stats, filter_fasta, make_workdir, next_step_name, require_databases, run_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(
    input: str,
    outfile: str,
    phylum: list[str],
    workdir: str | None = None,
    sourdb: str | None = None,
    sourdb_type: str = "gbk",
    kmer: str = "31",
    taxonomy: bool = False,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Classify each contig with sourmash and drop those classified outside ``phylum``.

    Each contig is classified with ``sourmash lca classify``; classified contigs whose lineage
    contains none of ``phylum`` are dropped, while unclassified (``nomatch``) contigs are kept.
    The sourmash CSV is copied next to ``input`` as ``<input stem>.sourmash-taxonomy.csv``.

    Args:
        input: Assembly FASTA to screen.
        outfile: Output FASTA of retained contigs.
        phylum: Taxon names (usually phyla) to keep.
        workdir: Working directory; a temporary one is created when None.
        sourdb: Path to a sourmash LCA database; when None the one for ``sourdb_type`` is used.
        sourdb_type: Installed database to use: ``gbk``, ``gtdbrep`` or ``gtdb``.
        kmer: k-mer size (as a string) used for ``sourmash compute``.
        taxonomy: Only report the taxonomic classifications and return without filtering.
        debug: Print per-contig details and keep the work directory.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.
    """
    workdir, custom_workdir = make_workdir(workdir, "taxpurge")
    try:
        if not sourdb:
            # --sourdb_type is restricted to these values by the "taxpurge" menu in _menu.py
            database = {"gbk": "sm_gbk", "gtdbrep": "sm_gtdbrep", "gtdb": "sm_gtdb"}[sourdb_type]
            sour_db = require_databases([database], hint="or pass --sourdb PATH")[0]
        else:
            sour_db = str(Path(sourdb).resolve())

        assembly_working = "assembly.fasta"
        shutil.copyfile(input, str(Path(workdir, assembly_working)))
        num_seqs, assembly_size = fasta_stats(str(Path(workdir, assembly_working)))
        logger.info(f"Assembly is {num_seqs:,} contigs and {assembly_size:,} bp")

        logger.info("Running SourMash to get taxonomy classification for each contig")
        sour_sketch = assembly_working + ".sig"
        sour_compute = ["sourmash", "compute", "-k", kmer, "--scaled=1000", "--singleton", assembly_working]
        run_cmd(sour_compute, debug, cwd=workdir)
        sour_classify = ["sourmash", "lca", "classify", "--db", sour_db, "--query", sour_sketch]
        # output csv: ID,status,superkingdom,phylum,class,order,family,genus,species,strain
        contig_taxonomy = {}
        unique_tax: set[str] = set()
        sourmash_tsv = str(Path(workdir, "sourmash.csv"))
        with open(sourmash_tsv, "w") as sour_out:
            for line in execute(sour_classify, workdir, debug):
                sour_out.write(line)
                if not line or line.startswith("\n") or line.startswith("ID") or line.count(",") < 9:
                    continue
                line = line.strip()
                cols = line.split(",")
                if cols[1].strip() == "found":
                    idx = 1
                    contig_taxonomy[cols[0]] = cols[idx + 1 :]
                    tax_clean = [x for x in cols[idx + 1 :] if x]
                    unique_tax.add(";".join(tax_clean))
                elif cols[1].strip() == "nomatch":
                    idx = cols.index("nomatch")
                    contig_taxonomy[cols[0]] = cols[idx + 1 :]
        logger.info("Found {:} taxonomic classifications for contigs:\n{:}".format(len(unique_tax), "\n".join(unique_tax)))
        _save_taxonomy_csv(sourmash_tsv, input)
        if taxonomy:  # --just-show-taxonomy only reports the classifications
            cleanup_workdir(workdir, debug, custom_workdir)
            return

        tax_to_drop = set()
        for k, v in contig_taxonomy.items():
            v = [x for x in v if x]  # remove empty items from list
            if debug:
                logger.info(f"{k}\t{v}")
            if v and not any(i in v for i in phylum):
                tax_to_drop.add(k)
        if debug:
            logger.info("Contigs dropped due to taxonomy: {:}".format(",".join(sorted(tax_to_drop))))

        logger.info(f"Dropping {len(tax_to_drop):,} contigs classified outside {', '.join(phylum)}")
        num_seqs, assembly_size = filter_fasta(str(Path(workdir, assembly_working)), outfile, lambda seq_id: seq_id not in tax_to_drop)
        logger.info(f"Taxonomy-purged assembly is {num_seqs:,} contigs and {assembly_size:,} bp")
        cleanup_workdir(workdir, debug, custom_workdir)

        if not pipe:
            next_out = next_step_name(outfile, ".covpurge.fasta")
            logger.info(f"Your next command might be:\nAAFTF covpurge -i {outfile} -o {next_out} -1 READ1 -2 READ2")
    finally:
        close_step_log()


def _save_taxonomy_csv(sourmash_tsv: str, input: str) -> None:
    """Copy the sourmash classification CSV next to the input as ``<input stem>.sourmash-taxonomy.csv``.

    Args:
        sourmash_tsv: The ``sourmash lca classify`` CSV in the work directory.
        input: The input assembly; the copy goes in its directory, named after it.
    """
    if not check_file(sourmash_tsv):
        return
    baseinput = Path(input).name
    if "." in baseinput:
        baseinput = baseinput.rsplit(".", 1)[0]
    shutil.copy(sourmash_tsv, str(Path(Path(input).parent, baseinput + ".sourmash-taxonomy.csv")))
