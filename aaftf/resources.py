"""Hard-coded URLs, accessions and database definitions used across AAFTF.

Holds the NCBI/EBI download locations for contaminant and vector databases, sourmash LCA
databases, NCBI FCS-adaptor release details, and the
``DATABASES`` table that ``AAFTF database`` downloads from.
"""

from typing import Any

__all__ = ["NCBI", "CONTAMINANT_ACCESSIONS", "DB_LINKS", "EUTILS", "SEQ_DBS", "FCSADAPTOR", "DATABASES"]


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
    "sourmash_gbk": [{"version": "2017.11.07", "filename": "genbank-k31.lca.json.gz", "url": "https://osf.io/4f8n3/download"}],
    # first in list is default, will fix someday to allow choosing the version
    # TODO: the GTDB entries below download sourmash .dna.zip signature collections but save them
    # under .lca.json.gz names, and taxpurge runs `sourmash lca classify --db`, which needs an LCA
    # database -- so --sourdb_type gtdb / gtdbrep likely fail. Point these at LCA databases (or
    # switch taxpurge to a classifier that takes zip collections). See TODO.md.
    "sourmash_gtdbrep": [
        {"version": "rs220", "filename": "gtdb-rs220-reps.k31.lca.json.gz", "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-reps-rs220-k31.dna.zip"},  # noqa: E501
    ],
    # first in list is default, will fix someday to allow choosing the version
    "sourmash_gtdb": [
        {"version": "rs220", "filename": "gtdb-rs220-k31.lca.json.gz", "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-rs220-k31.dna.zip"},  # noqa: E501
    ],
}

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
# read each one; `required` marks those `AAFTF pipeline_short` needs with its default settings (the rest are
# only for non-default options or optional steps); `executable` marks files that must be runnable
# after download.
DATABASES = {
    "phix": {"filename": CONTAMINANT_ACCESSIONS["phiX"][0].rsplit("/", 1)[-1], "url": CONTAMINANT_ACCESSIONS["phiX"][0], "used_by": "filter", "required": True},
    "univec": {"filename": DB_LINKS["UniVec"][0].rsplit("/", 1)[-1], "url": DB_LINKS["UniVec"][0], "used_by": "filter, vecscreen", "required": True},
    "euks": {"filename": DB_LINKS["CONTAM_EUKS"][0].rsplit("/", 1)[-1], "url": DB_LINKS["CONTAM_EUKS"][0], "used_by": "vecscreen", "required": True},
    "proks": {"filename": DB_LINKS["CONTAM_PROKS"][0].rsplit("/", 1)[-1], "url": DB_LINKS["CONTAM_PROKS"][0], "used_by": "vecscreen", "required": True},
    "mitodb": {"filename": DB_LINKS["MITO"][0].rsplit("/", 1)[-1], "url": DB_LINKS["MITO"][0], "used_by": "vecscreen", "required": True},
    "sm_gbk": {**{k: DB_LINKS["sourmash_gbk"][0][k] for k in ("filename", "url")}, "used_by": "taxpurge (--sourdb_type gbk)", "required": False},
    "sm_gtdbrep": {**{k: DB_LINKS["sourmash_gtdbrep"][0][k] for k in ("filename", "url")}, "used_by": "taxpurge (--sourdb_type gtdbrep)", "required": False},
    "sm_gtdb": {**{k: DB_LINKS["sourmash_gtdb"][0][k] for k in ("filename", "url")}, "used_by": "taxpurge (--sourdb_type gtdb)", "required": False},
    "fcs_script": {"filename": "run_fcsadaptor.sh", "url": FCSADAPTOR["EXEURL"] % FCSADAPTOR["VERSION"], "used_by": "fcs_screen", "required": False, "executable": True},
    "fcs_image": {
        "filename": FCSADAPTOR["SIFLOCAL"] % FCSADAPTOR["VERSION"],
        "url": "/".join([FCSADAPTOR["SIFURL"].rstrip("/"), FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"]]),
        "used_by": "fcs_screen (singularity)",
        "required": False,
    },
}
