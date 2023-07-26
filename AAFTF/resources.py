"""URLS and hardcoded resource links."""

NCBI = 'https://ftp.ncbi.nlm.nih.gov'

Contaminant_Accessions = {
    "phiX": [f'{NCBI}/genomes/all/GCF/000/819/615/' +
             'GCF_000819615.1_ViralProj14015/' +
             'GCF_000819615.1_ViralProj14015_genomic.fna.gz']
}


DB_Links = {
    'UniVec': [f'{NCBI}/pub/UniVec/UniVec'],
    'CONTAM_EUKS': [f'{NCBI}/pub/kitts/contam_in_euks.fa.gz'],
    'CONTAM_PROKS': [f'{NCBI}/pub/kitts/contam_in_prok.fa'],
    'MITO': [f'{NCBI}/refseq/release/mitochondrion/' +
             'mitochondrion.1.1.genomic.fna.gz'],
    'sourmash_gbk': [
        {'version': '2017.11.07',
         'filename': 'genbank-k31.lca.json.gz',
         'url': 'https://osf.io/4f8n3/download'}],
    # first in list is default, will fix someday to allow choosing the version
    'sourmash_gtdbrep': [
        {'version': 'rs214',
         'filename': 'tdb-rs214-reps.k31.lca.json.gz',
         'url': 'https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.lca.json.gz'},  # noqa: E501
        {'version': 'rs207',
         'filename': 'gtdb-rs207-genomic-reps.dna.k31.lca.json.gz',
         'url': 'https://osf.io/p9ezm/download'}],
    # first in list is default, will fix someday to allow choosing the version
    'sourmash_gtdb': [
        {'version': 'rs214',
         'filename': 'gtdb-rs214-k31.lca.json.gz',
         'url': 'https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.lca.json.gz'},  # noqa: E501
        {'version': 'rs207',
         'filename': 'gtdb-rs207.genomic.k31.lca.json.gz',
         'url': 'https://osf.io/tf3ah/download'}],
}

SeqDBs = {
    'nucleotide': [f'{NCBI}/entrez/eutils/efetch.fcgi?' +
                   'db=nucleotide&id=%s&rettype=fasta'],
    'nucleotide_ebi': 'https://www.ebi.ac.uk/ena/data/view/%s?display=fasta',
    'nucleotide_ncbi': [f'{NCBI}/entrez/eutils/efetch.fcgi?' +
                        'db=nucleotide&id=%s&rettype=fasta'],
}

Mitoseqs = {
    'COB1': ("atgagaattttaaaaagtcatcctttattaaaattagttaatagttatattattg" +
             "attcaccacaaccttctaatattagttatttatgaaattttggatctttattagc" +
             "tttatgtttagttatacaaattgtaactggtgttacattagctatgcactataca" +
             "cctaatgttgatttagcttttaattctgtagaacatattatgagagatgtaaata" +
             "atggttgattaataagatatttacatgctaatactgcttcagcattctttttctt" +
             "agttatatttacatataggtagaggattatattatggttcatataaatcacctag" +
             "aacttaacatgagctattgg"),
}

"""NCBI Foreign Contaminant Screen tool links"""
FCSADAPTOR = {
    'VERSION': '0.4.0',
    'SIF': 'fcs-adaptor.sif',
    'SIFLOCAL': 'fcs-adaptor.%s.sif',
    'SIFURL': 'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/',
    'EXEURL': 'https://raw.githubusercontent.com/ncbi/fcs/v%s/dist/run_fcsadaptor.sh',  # noqa: E501
}
