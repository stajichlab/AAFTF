
Contaminant_Accessions = {"phiX": 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna',
        }

DB_Links = {'UniVec': 'ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec',
            'CONTAM_EUKS': 'ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz',
            'CONTAM_PROKS': 'ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_prok.fa',
            'MITO': 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/mito.nt.gz',
        }

Applications = {'megablast': 'blastn',
                'trimmomatic': 'java -jar %s',
                }

SeqDBs = {'nucleotide': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta',
          'nucleotide_ebi': 'https://www.ebi.ac.uk/ena/data/view/%s?display=fasta',
          'nucleotide_ncbi': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta',
          'sourmash': 'https://osf.io/zskb9/download?version=1'
          }
