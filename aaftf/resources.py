"""Hard-coded URLs, accessions and database definitions used across AAFTF.

Holds the NCBI/EBI download locations for contaminant and vector databases, sourmash sketch
databases, NCBI FCS-adaptor release details, and the
``DATABASES`` table that ``AAFTF database`` downloads from.
"""

from typing import Any

__all__ = ["NCBI", "CONTAMINANT_ACCESSIONS", "DB_LINKS", "CONTAMDB", "CONTAMDB_BBTOOLS", "EUTILS", "SEQ_DBS", "FCSADAPTOR", "SOURMASH_DATABASES", "SOURMASH_LINEAGES", "TAXPURGE_DEFAULT_DATABASES", "DATABASES"]


# base URL of the NCBI FTP site
NCBI = "https://ftp.ncbi.nlm.nih.gov"

# genomes screened out of reads by ``filter`` (phiX spike-in)
CONTAMINANT_ACCESSIONS = {"phiX": [f"{NCBI}/genomes/all/GCF/000/819/615/" + "GCF_000819615.1_ViralProj14015/" + "GCF_000819615.1_ViralProj14015_genomic.fna.gz"]}

# NCBI's adaptor set for screening eukaryote assemblies (sequencing adaptors, primers and barcodes,
# 20-70 bp; no cloning vectors, so nothing in it can match fungal genomic DNA)
DB_LINKS: dict[str, list[Any]] = {"EUK_ADAPTORS": [f"{NCBI}/pub/kitts/adaptors_for_screening_euks.fa"]}

# contamination database shared by `filter` (reads) and `vecscreen` (contigs), built by
# `AAFTF database contamdb`: NCBI euk adaptors + PhiX + BBTools lambda.fa.gz and adapters.fa.
# `blastdb` is the path prefix of its BLAST database (contamdb.nhr, ...), next to the FASTA
CONTAMDB = {"fasta": "contamdb.fa", "blastdb": "contamdb"}
# BBTools resources/ files added to contamdb
CONTAMDB_BBTOOLS = ("lambda.fa.gz", "adapters.fa")

# sourmash databases (all k=31) for `taxpurge`, by `--sourdb` name. `AAFTF database` prepares each
# as a branchwater RocksDB index (folder `rocksdb`), which taxpurge searches from disk with much
# less memory: unpacked from a published archive (`rocksdb_url`) when there is one, else built
# locally from the sketch file (`filename` / `url`), which is then deleted (`--no-build` keeps the
# sketch file instead, and taxpurge then searches it). `lineages` is the lineage CSV (DATABASES
# abbreviation; None when taxpurge writes it itself); `scaled` is the sketches' resolution (one
# hash kept per `scaled` bp, so 10000 is 10x smaller but needs ~10x longer contigs)
_SM = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db"
SOURMASH_DATABASES: dict[str, dict[str, Any]] = {
    "gbk-fungi": {
        "lineages": "gbk-fungi-lineages",
        "scaled": 1000,
        "filename": "genbank-20250408-fungi-k31.zip",
        "url": f"{_SM}/genbank-2025.04/genbank-20250408-fungi-k31.zip",
        "rocksdb": "genbank-20250408-fungi-k31.rocksdb",
        "about": "GenBank 2025.04 fungi, all genomes",
    },
    "gbk-fungi-lite": {
        "lineages": "gbk-euks-lineages",
        "scaled": 10000,
        "filename": "genbank-euks-2025.01-fungi.k31.sig.zip",
        "url": f"{_SM}/genbank-euks-2025.01/fungi.k31.sig.zip",
        "rocksdb": "genbank-euks-2025.01-fungi.k31.rocksdb",
        "about": "GenBank 2025.01 eukaryote collection, fungi",
    },
    "gtdb-reps": {
        "lineages": "gtdb-reps-lineages",
        "scaled": 1000,
        "rocksdb": "gtdb-rs226-reps.k31.rocksdb",
        "rocksdb_url": f"{_SM}/gtdb-rs226/gtdb-rs226-reps.k31.rocksdb.tar.gz",
        "about": "GTDB rs226 representative genomes, bacteria and archaea",
    },
    "gtdb-reps-lite": {
        "lineages": "gtdb-reps-lineages",
        "scaled": 10000,
        "rocksdb": "gtdb-rs226-reps.k31-sc10k.sig.rocksdb",
        "rocksdb_url": f"{_SM}/gtdb-rs226/gtdb-rs226-reps.k31-sc10k.sig.rocksdb.tar.gz",
        "about": "GTDB rs226 representative genomes, bacteria and archaea",
    },
    "human-hg38": {"lineages": None, "scaled": 1000, "filename": "hg38-entire.sig.zip", "url": f"{_SM}/hg38/hg38-entire.sig.zip", "rocksdb": "hg38-entire.rocksdb", "about": "human reference genome hg38"},
}
SOURMASH_LINEAGES = {
    "gbk-fungi-lineages": {"filename": "genbank-20250408-lineages.fungi.csv", "url": f"{_SM}/genbank-2025.04/lineages.fungi.csv"},
    "gbk-euks-lineages": {"filename": "genbank-euks-2025.01-eukaryotes.lineages.csv", "url": f"{_SM}/genbank-euks-2025.01/eukaryotes.lineages.csv"},
    "gtdb-reps-lineages": {"filename": "gtdb-rs226-reps.lineages.csv", "url": f"{_SM}/gtdb-rs226/gtdb-rs226-reps.lineages.csv"},
}


def _sourmash_database_entries() -> dict[str, dict[str, Any]]:
    """Return the DATABASES entries of the sourmash databases, each followed by its lineage CSV.

    A database with a published index (``rocksdb_url``) is that index folder (``index``); one
    without is its sketch file plus the index folder ``AAFTF database`` builds from it
    (``build_index``).
    """
    entries: dict[str, dict[str, Any]] = {}
    for name, entry in SOURMASH_DATABASES.items():
        common = {"used_by": f"taxpurge (--sourdb {name})", "required": False, "with": [entry["lineages"]] if entry["lineages"] else []}
        if entry.get("rocksdb_url"):
            entries[name] = {"filename": entry["rocksdb"], "url": entry["rocksdb_url"], "index": True, **common}
        else:
            entries[name] = {"filename": entry["filename"], "url": entry["url"], "build_index": entry["rocksdb"], "scaled": entry["scaled"], **common}
        lineage = entry["lineages"]
        if lineage and lineage not in entries:
            users = ", ".join(db for db, e in SOURMASH_DATABASES.items() if e["lineages"] == lineage)
            entries[lineage] = {**SOURMASH_LINEAGES[lineage], "used_by": f"taxpurge (lineages of {users})", "required": False}
    return entries


# `taxpurge --sourdb` default (also downloaded by `AAFTF database taxpurge`)
TAXPURGE_DEFAULT_DATABASES = ["gbk-fungi-lite", "gtdb-reps-lite", "human-hg38"]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# URL templates for fetching one sequence by accession (``% accession``)
SEQ_DBS = {
    "nucleotide": f"{EUTILS}/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta",
    "nucleotide_ebi": "https://www.ebi.ac.uk/ena/data/view/%s?display=fasta",
    "nucleotide_ncbi": f"{EUTILS}/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta",
}

# NCBI Foreign Contaminant Screen (FCS-adaptor) release, container image and script links

FCSADAPTOR = {
    "VERSION": "0.5.5",
    "SIF": "fcs-adaptor.sif",
    "SIFLOCAL": "fcs-adaptor.%s.sif",
    "SIFURL": "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/",
    "EXEURL": "https://raw.githubusercontent.com/ncbi/fcs/v%s/dist/run_fcsadaptor.sh",  # noqa: E501
    "DOCKERIMAGE": "ncbi/fcs-adaptor:%s",
}

# Databases `AAFTF database` can fetch, keyed by abbreviation. `used_by` lists the subcommands that
# read each one; `with` lists companion databases downloaded along with it; `index` marks a RocksDB
# index folder (downloaded as a .tar.gz archive and unpacked) and `build_index` the index folder
# built from a sketch file (see SOURMASH_DATABASES); `required` marks those `AAFTF pipeline_short` needs with its default settings (the rest are
# only for non-default options or optional steps); `executable` marks files that must be runnable
# after download; `built_from` marks a database `AAFTF database` builds locally from the listed
# databases (no `url`).
DATABASES = {
    "phix": {"filename": CONTAMINANT_ACCESSIONS["phiX"][0].rsplit("/", 1)[-1], "url": CONTAMINANT_ACCESSIONS["phiX"][0], "used_by": "contamdb (source)", "required": True},
    "euk_adaptors": {"filename": DB_LINKS["EUK_ADAPTORS"][0].rsplit("/", 1)[-1], "url": DB_LINKS["EUK_ADAPTORS"][0], "used_by": "contamdb (source)", "required": True},
    "contamdb": {"filename": CONTAMDB["fasta"], "built_from": ["euk_adaptors", "phix"], "used_by": "filter, vecscreen", "required": True},
    **_sourmash_database_entries(),
    "fcs_script": {"filename": "run_fcsadaptor.sh", "url": FCSADAPTOR["EXEURL"] % FCSADAPTOR["VERSION"], "used_by": "fcs_screen", "required": False, "executable": True},
    "fcs_image": {
        "filename": FCSADAPTOR["SIFLOCAL"] % FCSADAPTOR["VERSION"],
        "url": "/".join([FCSADAPTOR["SIFURL"].rstrip("/"), FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"]]),
        "used_by": "fcs_screen (singularity)",
        "required": False,
    },
}
