import sys, os, shutil

from subprocess import call, Popen, PIPE, STDOUT
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from AAFTF.utility import execute
from AAFTF.utility import calcN50
from AAFTF.utility import countfasta
from AAFTF.resources import SeqDBs

#logging

import logging
logger = logging.getLogger('AAFTF')


def run(parser,args):

    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    
    #find reads
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        for file in os.listdir(args.tmpdir):
            if '_cleaned' in file and file.endswith('.fq.gz'):
                if '_1' in file:
                    forReads = os.path.abspath(os.path.join(args.tmpdir, file))
                if '_2' in file:
                    revReads = os.path.abspath(os.path.join(args.tmpdir, file))
    if not forReads:
        print('Unable to located FASTQ raw reads')
        sys.exit(1)
        
    # hard coded tmpfile
    assembly_working  = 'assembly.fasta'
    megablast_working = 'megablast.out'
    blobBAM           = "remapped.bam"
    shutil.copyfile(args.input, os.path.join(args.tmpdir,assembly_working))
    logger.info('Loaded assembly containing {:,} contigs'.format(countfasta(os.path.join(args.tmpdir,assembly_working))))
    DEVNULL = open(os.devnull, 'w')
    #check if BAM present, if so skip 
    if not os.path.isfile(os.path.join(args.tmpdir, blobBAM)):  
        # index
        bwa_index  = ['bwa','index', os.path.basename(assembly_working)]     
        if not os.path.isfile(os.path.join(args.tmpdir,os.path.basename(assembly_working)+".amb")):
            logger.info('Building BWA index')
            logger.info('CMD: {:}'.format(' '.join(bwa_index)))
            subprocess.run(bwa_index, cwd=args.tmpdir, stderr=DEVNULL)
        else:
            logger.info("BWA index already exists for %s"%(os.path.join(args.tmpdir,os.path.basename(assembly_working))))
        #mapped reads to assembly using BWA
        bwa_cmd = ['bwa','mem',
                   '-t', str(args.cpus),
                    os.path.basename(assembly_working), # assembly index base
                    forReads]    
        if revReads:
            bwa_cmd.append(revReads)
    
        #run BWA and pipe to samtools sort
        logger.info('Aligning reads to assembly with BWA')
        logger.info('CMD: {:}'.format(' '.join(bwa_cmd)))
        p1 = subprocess.Popen(bwa_cmd, cwd=args.tmpdir,
                              stdout=subprocess.PIPE, stderr=DEVNULL)
        p2 = subprocess.Popen(['samtools', 'sort', '--threads', str(args.cpus),
                               '-o', blobBAM, '-'], cwd=args.tmpdir,
                              stdout=subprocess.PIPE, stderr=DEVNULL,
                              stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()

        if os.path.exists(blobBAM):
            subprocess.run(['samtools', 'index', blobBAM], cwd=args.tmpdir)
    
    #now calculate coverage
    logger.info('Calculating read coverage per contig')
    FastaBed = os.path.join(args.tmpdir, 'assembly.bed')
    lengths = []
    with open(FastaBed, 'w') as bedout:
        with open(os.path.join(args.tmpdir, assembly_working), 'rU') as SeqIn:
            for record in SeqIO.parse(SeqIn, 'fasta'):
                bedout.write('{:}\t{:}\t{:}\n'.format(record.id, 0, len(record.seq)))
                lengths.append(len(record.seq))
    
    N50 = calcN50(lengths)
    Coverage = {}
    cov_cmd = ['samtools', 'bedcov', os.path.basename(FastaBed), blobBAM] 
    logger.info('CMD: {:}'.format(' '.join(cov_cmd)))
    for line in execute(cov_cmd, args.tmpdir):
        if not line or line.startswith('\n') or line.count('\t') < 3:
            continue
        line = line.strip()
        cols = line.split('\t')
        cov = int(cols[3]) / float(cols[2])
        Coverage[cols[0]] = (int(cols[2]), cov)
    
    #get average coverage of N50 contigs
    n50Cov = []
    for k,v in Coverage.items():
        if args.debug:
            print('{:}; Len: {:}; Cov: {:.2f}'.format(k, v[0], v[1]))
        if v[0] >= N50:
            n50Cov.append(v[1])
    n50AvgCov = sum(n50Cov) / len(n50Cov)
    minpct = args.mincovpct / 100
    min_coverage = int(n50AvgCov * minpct)  #should we make this a variable? 5% was something arbitrary
    logger.info('Average coverage for N50 contigs is {:}X'.format(int(n50AvgCov)))
    #Start list of contigs to drop
    Contigs2Drop = []
    for k,v in Coverage.items():
        if v[1] <= min_coverage:
            Contigs2Drop.append(k)
    logger.info('Found {:,} contigs with coverage less than {:}X ({:}%)'.format(len(Contigs2Drop), 
    			min_coverage, args.mincovpct))
            
    #now filter for taxonomy with sourmash lca classify
    # sourmash compute -k 31 --scaled=1000 --singleton NRRL_66455_July2018.polished_assembly.fasta 
    # sourmash lca classify --db genbank-k31.lca.json.gz --query NRRL_66455_July2018.polished_assembly.fasta.sig -o output.csv
    logger.info('Running SourMash to get taxonomy classification for each contig')
    sour_sketch = os.path.basename(assembly_working)+'.sig'
    sour_compute = ['sourmash', 'compute', '-k', '31', '--scaled=1000',
                   '--singleton', assembly_working]
    logger.info('CMD: {:}'.format(' '.join(sour_compute)))
    subprocess.run(sour_compute, cwd=args.tmpdir, stderr=DEVNULL)
    
    sour_classify = ['sourmash', 'lca', 'classify', '--db', os.path.abspath(args.sourdb),
                     '--query', sour_sketch]
    logger.info('CMD: {:}'.format(' '.join(sour_classify)))
    # output csv: ID,status,superkingdom,phylum,class,order,family,genus,species,strain
    Taxonomy = {}
    UniqueTax = []
    for line in execute(sour_classify, args.tmpdir):
        if not line or line.startswith('\n') or line.startswith('ID') or line.count(',') < 9:
            continue
        line = line.strip()
        cols = line.split(',')
        if 'found' in cols:
            idx = cols.index('found')
            Taxonomy[cols[0]] = cols[idx+1:]
            taxClean = [x for x in cols[idx+1:] if x]
            UniqueTax.append('{:}'.format(';'.join(taxClean)))
        elif 'nomatch' in cols:
            idx = cols.index('nomatch')
            Taxonomy[cols[0]] = cols[idx+1:]
    UniqueTax = set(UniqueTax)
    logger.info('Found {:} taxonomic classifications for contigs:\n{:}'.format(len(UniqueTax), '\n'.join(UniqueTax)))
    Tax2Drop = []
    for k,v in Taxonomy.items():
        v = [x for x in v if x] #remove empty items from list
        if args.debug:
            print('{:}\t{:}'.format(k, v))
        if len(v) > 0:
            if not any(i in v for i in args.phylum):
                Tax2Drop.append(k)
    logger.info('{:,} contigs are contamination based on taxonomy'.format(len(Tax2Drop)))
    if args.debug:
        print('Contigs dropped due to coverage: {:}'.format(','.join(Contigs2Drop)))
        print('Contigs dropped due to taxonomy: {:}'.format(','.join(Tax2Drop)))
    DropFinal = Contigs2Drop + Tax2Drop
    DropFinal = set(DropFinal)
    logger.info('Dropping {:} contigs'.format(len(DropFinal)))
    with open(args.out, 'w') as outfile:
        with open(os.path.join(args.tmpdir,assembly_working), 'rU') as seqin:
            for record in SeqIO.parse(seqin, 'fasta'):
                if not record.id in DropFinal:
                    SeqIO.write(record, outfile, 'fasta')
    logger.info('Cleaned assembly containing {:,} contigs: {:}'.format(countfasta(args.out), args.out))
