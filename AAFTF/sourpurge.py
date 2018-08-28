import sys, os, shutil

from subprocess import call, Popen, PIPE, STDOUT
import subprocess

from Bio import SeqIO
from AAFTF.utility import execute
from AAFTF.utility import calcN50
from AAFTF.utility import fastastats
from AAFTF.resources import DB_Links

# logging - we may need to think about whether this has 
# separate name for the different runfolder

import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):

    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    if args.prefix:
        prefix = args.prefix
    elif args.left:
        prefix = os.path.basename(os.path.splitext(args.left)[0])
    elif args.outfile:
        prefix = os.path.basename(os.path.splitext(args.outfile)[0])

    if not prefix.endswith('_'):
        prefix += "_"

    #find reads
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        for file in os.listdir(args.workdir):
            if '_cleaned' in file and file.endswith('q.gz') and file.startswith(prefix):
                if '_1.fastq' in file:
                    forReads = os.path.abspath(os.path.join(args.workdir, file))
                if '_2.fastq' in file:
                    revReads = os.path.abspath(os.path.join(args.workdir, file))
    if not forReads:
        print('Unable to located FASTQ raw reads')
        sys.exit(1)
    
    #parse database locations
    if not args.sourdb:
        try:
            DB = os.environ["AAFTF_DB"]
        except KeyError:
            if args.AAFTF_DB:
                SOUR = os.path.join(args.AAFTF_DB, 'genbank-k31.lca.json.gz')
            else:
                logger.error("$AAFTF_DB/genbank-k31.lca.json.gz not found, pass --sourdb")
                sys.exit(1)
        SOUR = os.path.join(DB, 'genbank-k31.lca.json.gz')
        if not os.path.isfile(SOUR):
            logger.error("{:} sourmash database not found".format(SOUR))
            # should we prompt it to download 
            exit
    else:
        SOUR = os.path.abspath(args.sourdb)
                    
    # hard coded tmpfile
    assembly_working  = prefix + 'assembly.fasta'
    megablast_working = prefix + 'megablast.out'
    blobBAM           = prefix + 'remapped.bam'
    shutil.copyfile(args.input, os.path.join(args.workdir,assembly_working))
    numSeqs, assemblySize = fastastats(os.path.join(args.workdir,assembly_working))
    logger.info('Assembly is {:,} contigs and {:,} bp'.format(numSeqs, assemblySize))
    DEVNULL = open(os.devnull, 'w')

    #now filter for taxonomy with sourmash lca classify
    # sourmash compute -k 31 --scaled=1000 --singleton NRRL_66455_July2018.polished_assembly.fasta 
    # sourmash lca classify --db genbank-k31.lca.json.gz --query NRRL_66455_July2018.polished_assembly.fasta.sig -o output.csv
    logger.info('Running SourMash to get taxonomy classification for each contig')
    sour_sketch = os.path.basename(assembly_working)+'.sig'
    sour_compute = ['sourmash', 'compute', '-k', '31', '--scaled=1000',
                   '--singleton', assembly_working]
    logger.info('CMD: {:}'.format(' '.join(sour_compute)))
    subprocess.run(sour_compute, cwd=args.workdir, stderr=DEVNULL)
    
    sour_classify = ['sourmash', 'lca', 'classify', '--db', SOUR,
                     '--query', sour_sketch]
    logger.info('CMD: {:}'.format(' '.join(sour_classify)))
    # output csv: ID,status,superkingdom,phylum,class,order,family,genus,species,strain
    Taxonomy = {}
    UniqueTax = []
    sourmashTSV = os.path.join(args.workdir, prefix + 'sourmash.tsv')
    with open(sourmashTSV, 'w') as sour_out:
        for line in execute(sour_classify, args.workdir):
            sour_out.write(line)
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
    if args.taxonomy:
        sys.exit(1)
    Tax2Drop = []
    for k,v in Taxonomy.items():
        v = [x for x in v if x] #remove empty items from list
        if args.debug:
            print('{:}\t{:}'.format(k, v))
        if len(v) > 0:
            if not any(i in v for i in args.phylum):
                Tax2Drop.append(k)
    
    #drop contigs from taxonomy before calculating coverage
    logger.info('Dropping {:} contigs from taxonomy screen'.format(len(Tax2Drop)))
    sourTax = os.path.join(args.workdir, prefix+'sourmashed-tax-screen.fasta')
    with open(sourTax, 'w') as outfile:
        with open(os.path.join(args.workdir,assembly_working), 'rU') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                if not record.id in Tax2Drop:
                    SeqIO.write(record, outfile, 'fasta')
        
    #check if BAM present, if so skip 
    if not os.path.isfile(os.path.join(args.workdir, blobBAM)):  
        # index
        bwa_index  = ['bwa','index', os.path.basename(sourTax)]     
        logger.info('Building BWA index')
        logger.info('CMD: {:}'.format(' '.join(bwa_index)))
        subprocess.run(bwa_index, cwd=args.workdir, stderr=DEVNULL)
        #mapped reads to assembly using BWA
        bwa_cmd = ['bwa','mem',
                   '-t', str(args.cpus),
                    os.path.basename(sourTax), # assembly index base
                    forReads]    
        if revReads:
            bwa_cmd.append(revReads)
    
        #run BWA and pipe to samtools sort
        logger.info('Aligning reads to assembly with BWA')
        logger.info('CMD: {:}'.format(' '.join(bwa_cmd)))
        p1 = subprocess.Popen(bwa_cmd, cwd=args.workdir,
                              stdout=subprocess.PIPE, stderr=DEVNULL)
        p2 = subprocess.Popen(['samtools', 'sort', '--threads', str(args.cpus),
                               '-o', blobBAM, '-'], cwd=args.workdir,
                              stdout=subprocess.PIPE, stderr=DEVNULL,
                              stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()

        subprocess.run(['samtools', 'index', blobBAM], cwd=args.workdir)


    #now calculate coverage
    logger.info('Calculating read coverage per contig')
    FastaBed = os.path.join(args.workdir, prefix+'assembly.bed')
    lengths = []
    with open(FastaBed, 'w') as bedout:
        with open(sourTax, 'rU') as SeqIn:
            for record in SeqIO.parse(SeqIn, 'fasta'):
                bedout.write('{:}\t{:}\t{:}\n'.format(record.id, 0, len(record.seq)))
                lengths.append(len(record.seq))
    
    N50 = calcN50(lengths)
    Coverage = {}
    coverageBed = os.path.join(args.workdir, prefix+'coverage.bed')
    cov_cmd = ['samtools', 'bedcov', os.path.basename(FastaBed), blobBAM] 
    logger.info('CMD: {:}'.format(' '.join(cov_cmd)))
    with open(coverageBed, 'w') as bed_out:
        for line in execute(cov_cmd, args.workdir):
            bed_out.write(line)
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
    min_coverage = float(n50AvgCov * minpct)  #should we make this a variable? 5% was something arbitrary
    logger.info('Average coverage for N50 contigs is {:}X'.format(int(n50AvgCov)))
    #Start list of contigs to drop
    Contigs2Drop = []
    for k,v in Coverage.items():
        if v[1] <= min_coverage:
            Contigs2Drop.append(k)
    logger.info('Found {:,} contigs with coverage less than {:.2f}X ({:}%)'.format(len(Contigs2Drop), 
                min_coverage, args.mincovpct))
            
    if args.debug:
        print('Contigs dropped due to coverage: {:}'.format(','.join(Contigs2Drop)))
        print('Contigs dropped due to taxonomy: {:}'.format(','.join(Tax2Drop)))
    DropFinal = Contigs2Drop + Tax2Drop
    DropFinal = set(DropFinal)
    logger.info('Dropping {:} total contigs based on taxonomy and coverage'.format(len(DropFinal)))
    with open(args.outfile, 'w') as outfile:
        with open(sourTax, 'rU') as seqin:
            for record in SeqIO.parse(seqin, 'fasta'):
                if not record.id in DropFinal:
                    SeqIO.write(record, outfile, 'fasta')
                    
    numSeqs, assemblySize = fastastats(args.outfile)
    logger.info('Cleaned assembly is {:,} contigs and {:,} bp'.format(numSeqs, assemblySize))
    if '_' in args.outfile:
        nextOut = args.outfile.split('_')[0]+'.rmdup.fasta'
    elif '.' in args.out:
        nextOut = args.outfile.split('.')[0]+'.rmdup.fasta'
    else:
        nextOut = args.outfile+'.rmdup.fasta'
    
    logger.info('Your next command might be:\n\tAAFTF rmdup -i {:} -o {:}\n'.format(
                args.outfile, nextOut))
