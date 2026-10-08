"""Hard-coded URLs, accessions and database definitions used across AAFTF.

Holds the NCBI/EBI download locations for contaminant and vector databases, sourmash sketch
databases, NCBI FCS-adaptor release details, and the
``DATABASES`` table that ``AAFTF database`` downloads from.
"""

from typing import Any

__all__ = ["NCBI", "CONTAMINANT_ACCESSIONS", "DB_LINKS", "EUTILS", "SEQ_DBS", "FCSADAPTOR", "SOURMASH_DATABASES", "SOURMASH_LINEAGES", "TAXPURGE_DEFAULT_DATABASES", "DATABASES"]


# base URL of the NCBI FTP site
NCBI = "https://ftp.ncbi.nlm.nih.gov"

# genomes screened out of reads by ``filter`` (phiX spike-in)
CONTAMINANT_ACCESSIONS = {"phiX": [f"{NCBI}/genomes/all/GCF/000/819/615/" + "GCF_000819615.1_ViralProj14015/" + "GCF_000819615.1_ViralProj14015_genomic.fna.gz"]}

# vector/contaminant/organelle FASTA URLs, and sourmash databases as {version, filename, url} dicts
DB_LINKS: dict[str, list[Any]] = {
    "UniVec": [f"{NCBI}/pub/UniVec/UniVec"],
    "CONTAM_EUKS": [f"{NCBI}/pub/kitts/contam_in_euks.fa.gz"],
    "CONTAM_PROKS": [f"{NCBI}/pub/kitts/contam_in_prok.fa"],
    "MITO": [f"{NCBI}/refseq/release/mitochondrion/" + "mitochondrion.1.1.genomic.fna.gz"],
}

# sourmash sketch collections (all k=31) for `taxpurge`, by `--sourdb` name: the sketch file, its
# lineage CSV (DATABASES abbreviation; None when taxpurge writes it itself) and the sketches'
# `scaled` (one hash kept per `scaled` bp, so 10000 is 10x smaller but needs ~10x longer contigs)
_SM = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db"
SOURMASH_DATABASES: dict[str, dict[str, Any]] = {
    "gbk-fungi": {"lineages": "gbk-fungi-lineages", "scaled": 1000, "filename": "genbank-20250408-fungi-k31.zip", "url": f"{_SM}/genbank-2025.04/genbank-20250408-fungi-k31.zip", "about": "GenBank 2025.04 fungi, all genomes (5.3 GB)"},
    "gbk-fungi-lite": {"lineages": "gbk-euks-lineages", "scaled": 10000, "filename": "genbank-euks-2025.01-fungi.k31.sig.zip", "url": f"{_SM}/genbank-euks-2025.01/fungi.k31.sig.zip", "about": "GenBank 2025.01 eukaryote collection, fungi (170 MB)"},
    "gtdb-reps": {"lineages": "gtdb-reps-lineages", "scaled": 1000, "filename": "gtdb-rs226-reps.k31.sig.zip", "url": f"{_SM}/gtdb-rs226/gtdb-rs226-reps.k31.sig.zip", "about": "GTDB rs226 representative genomes, bacteria and archaea (3.9 GB)"},
    "gtdb-reps-lite": {"lineages": "gtdb-reps-lineages", "scaled": 10000, "filename": "gtdb-rs226-reps.k31-sc10k.sig.zip", "url": f"{_SM}/gtdb-rs226/gtdb-rs226-reps.k31-sc10k.sig.zip", "about": "GTDB rs226 representative genomes, bacteria and archaea (447 MB)"},
    "human-hg38": {"lineages": None, "scaled": 1000, "filename": "hg38-entire.sig.zip", "url": f"{_SM}/hg38/hg38-entire.sig.zip", "about": "human reference genome hg38 (55 MB)"},
}
SOURMASH_LINEAGES = {
    "gbk-fungi-lineages": {"filename": "genbank-20250408-lineages.fungi.csv", "url": f"{_SM}/genbank-2025.04/lineages.fungi.csv"},
    "gbk-euks-lineages": {"filename": "genbank-euks-2025.01-eukaryotes.lineages.csv", "url": f"{_SM}/genbank-euks-2025.01/eukaryotes.lineages.csv"},
    "gtdb-reps-lineages": {"filename": "gtdb-rs226-reps.lineages.csv", "url": f"{_SM}/gtdb-rs226/gtdb-rs226-reps.lineages.csv"},
}
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
# read each one; `with` lists companion databases downloaded along with it; `required` marks those `AAFTF pipeline_short` needs with its default settings (the rest are
# only for non-default options or optional steps); `executable` marks files that must be runnable
# after download.
DATABASES = {
    "phix": {"filename": CONTAMINANT_ACCESSIONS["phiX"][0].rsplit("/", 1)[-1], "url": CONTAMINANT_ACCESSIONS["phiX"][0], "used_by": "filter", "required": True},
    "univec": {"filename": DB_LINKS["UniVec"][0].rsplit("/", 1)[-1], "url": DB_LINKS["UniVec"][0], "used_by": "filter, vecscreen", "required": True},
    "euks": {"filename": DB_LINKS["CONTAM_EUKS"][0].rsplit("/", 1)[-1], "url": DB_LINKS["CONTAM_EUKS"][0], "used_by": "vecscreen", "required": True},
    "proks": {"filename": DB_LINKS["CONTAM_PROKS"][0].rsplit("/", 1)[-1], "url": DB_LINKS["CONTAM_PROKS"][0], "used_by": "vecscreen", "required": True},
    "mitodb": {"filename": DB_LINKS["MITO"][0].rsplit("/", 1)[-1], "url": DB_LINKS["MITO"][0], "used_by": "vecscreen", "required": True},
    **{name: {"filename": entry["filename"], "url": entry["url"], "used_by": f"taxpurge (--sourdb {name})", "required": False, "with": [entry["lineages"]] if entry["lineages"] else []} for name, entry in SOURMASH_DATABASES.items()},
    **{name: {**entry, "used_by": "taxpurge (lineages of " + ", ".join(db for db, e in SOURMASH_DATABASES.items() if e["lineages"] == name) + ")", "required": False} for name, entry in SOURMASH_LINEAGES.items()},
    "fcs_script": {"filename": "run_fcsadaptor.sh", "url": FCSADAPTOR["EXEURL"] % FCSADAPTOR["VERSION"], "used_by": "fcs_screen", "required": False, "executable": True},
    "fcs_image": {
        "filename": FCSADAPTOR["SIFLOCAL"] % FCSADAPTOR["VERSION"],
        "url": "/".join([FCSADAPTOR["SIFURL"].rstrip("/"), FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"]]),
        "used_by": "fcs_screen (singularity)",
        "required": False,
    },
}
