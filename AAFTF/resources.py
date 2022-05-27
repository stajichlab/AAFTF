
Contaminant_Accessions = {"phiX": 'https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna',
                          }

DB_Links = {'UniVec': 'https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec',
            'CONTAM_EUKS': 'https://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz',
            'CONTAM_PROKS': 'https://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_prok.fa',
            'MITO': 'https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz',
            # 'sourmash': 'https://osf.io/9xdg2/download?version=1'
            # can store multiple versions here
            'sourmash_gbk': [{'version': 2018,
                              'filename': 'genbank-k31.lca.json.gz',
                              'url': 'https://osf.io/p9ezm/download'}],
            'sourmash_gtdb_rep': [{'version': 'rs207',
                                   'filename': 'gtdb-rs207-genomic-reps.dna.k31.lca.json.gz',
                                   'url':      'https://osf.io/p9ezm/download'}],
            'sourmash_gtdb': [{'version': 'rs207',
                               'filename': 'gtdb-rs207.genomic.k31.lca.json.gz',
                               'url':      'https://osf.io/tf3ah/download'}]
            }


SeqDBs = {'nucleotide': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta',
          'nucleotide_ebi': 'https://www.ebi.ac.uk/ena/data/view/%s?display=fasta',
          'nucleotide_ncbi': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta',
          }
