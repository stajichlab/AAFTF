"""URLS and hardcoded resource links."""

NCBI = 'https://ftp.ncbi.nlm.nih.gov'

Contaminant_Accessions = {
    "phiX": (NCBI + '/genomes/Viruses/' +
             'enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna'),
}


DB_Links = {
    'UniVec': [NCBI+'/pub/UniVec/UniVec'],
    'CONTAM_EUKS': [NCBI + '/pub/kitts/contam_in_euks.fa.gz'],
    'CONTAM_PROKS': [NCBI + '/pub/kitts/contam_in_prok.fa'],
    'MITO': [NCBI + '/refseq/release/mitochondrion/' +
             'mitochondrion.1.1.genomic.fna.gz',
             NCBI + '/refseq/release/mitochondrion/' +
             'mitochondrion.2.1.genomic.fna.gz'],
    # 'sourmash': 'https://osf.io/9xdg2/download?version=1'
    # can store multiple versions here
    'sourmash_gbk': [{
        'version': 2018,
        'filename': 'genbank-k31.lca.json.gz',
        'url': 'https://osf.io/p9ezm/download'}],
    'sourmash_gtdbrep': [{
        'version': 'rs207',
        'filename': 'gtdb-rs207-genomic-reps.dna.k31.lca.json.gz',
        'url':      'https://osf.io/p9ezm/download'}],
    'sourmash_gtdb': [{
        'version': 'rs207',
        'filename': 'gtdb-rs207.genomic.k31.lca.json.gz',
        'url':      'https://osf.io/tf3ah/download'}]
}

SeqDBs = {
    'nucleotide': (NCBI + '/entrez/eutils/efetch.fcgi?' +
                   'db=nucleotide&id=%s&rettype=fasta'),
    'nucleotide_ebi': 'https://www.ebi.ac.uk/ena/data/view/%s?display=fasta',
    'nucleotide_ncbi': (NCBI + '/entrez/eutils/efetch.fcgi?' +
                        'db=nucleotide&id=%s&rettype=fasta'),
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
    'EXEURL': 'https://raw.githubusercontent.com/ncbi/fcs/v%s/dist/run_fcsadaptor.sh',
}
