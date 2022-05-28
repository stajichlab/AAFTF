import sys, os, shutil
from subprocess import call, Popen, PIPE, STDOUT
import subprocess
import uuid
from Bio import SeqIO
from AAFTF.utility import execute
from AAFTF.utility import calcN50
from AAFTF.utility import fastastats
from AAFTF.resources import DB_Links
from AAFTF.utility import status
from AAFTF.utility import printCMD
from AAFTF.utility import SafeRemove
from AAFTF.utility import checkfile
from AAFTF.utility import download

# logging - we may need to think about whether this has
# separate name for the different runfolder
def run(parser,args):

    if not args.workdir:
        args.workdir = 'aaftf-sourpurge_'+str(uuid.uuid4())[:8]
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    bamthreads = 4
    if args.cpus < 4:
        bamthreads = 1

    #find reads
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        status('Unable to located FASTQ raw reads, low coverage will be skipped. Provide -l,--left or -r,--right to enable low coverage filtering.')
#        sys.exit(1)

    #parse database locations
    if not args.sourdb:
#        'genbank-k31.lca.json.gz'
        dbfle="genbank-k31.lca.json.gz" # old default

        if args.sourdb_type.lower() == "gtdb":
            dbfile = DB_Links['sourmash_gtdb'][0]['filename']
        elif args.sourdb_type.lower() == "gtdbrep" or args.sourdb_type.lower() == "gtdb_rep":
            dbfile = DB_Links['sourmash_gtdbrep'][0]['filename']
        elif args.sourdb_type.lower() == "gbk" or args.sourdb_type.lower() = "genbank":
            dbfile = DB_Links['sourmash_gbk'][0]['filename']
        else:
            status("Unknown sourdb_type value {:} use one of {}".format(args.sourdb_type, ['gtdb','gtdbrep','gbk']))
            sys.exit(1)

        try:
            DB = os.environ["AAFTF_DB"]
        except KeyError:
            if args.AAFTF_DB:
                SOUR = os.path.join(args.AAFTF_DB, dbfile)
            else:
                status("$AAFTF_DB/{} not found, pass --sourdb".format(dbfile))
                sys.exit(1)
        SOUR = os.path.join(DB, dbfile)
        if not os.path.isfile(SOUR):
            status("{:} sourmash database not found, download and rename to {}/{}".format(SOUR,args.AAFTF_DB,dbfile))
            sys.exit(1)
    else:
        SOUR = os.path.abspath(args.sourdb)

    # hard coded tmpfile
    assembly_working  = 'assembly.fasta'
    blobBAM           = 'remapped.bam'
    shutil.copyfile(args.input, os.path.join(args.workdir,assembly_working))
    numSeqs, assemblySize = fastastats(os.path.join(args.workdir,
                                                    assembly_working))
    status('Assembly is {:,} contigs and {:,} bp'.format(numSeqs,
                                                              assemblySize))
    DEVNULL = open(os.devnull, 'w')

    #now filter for taxonomy with sourmash lca classify
    status('Running SourMash to get taxonomy classification for each contig')
    sour_sketch = os.path.basename(assembly_working)+'.sig'

    sour_compute = ['sourmash', 'compute', '-k', '31', '--scaled=1000',
                   '--singleton', assembly_working]
    printCMD(sour_compute)
    subprocess.run(sour_compute, cwd=args.workdir, stderr=DEVNULL)
    sour_classify = ['sourmash', 'lca', 'classify', '--db', SOUR,'--query', sour_sketch]
    printCMD(sour_classify)
    # output csv: ID,status,superkingdom,phylum,class,order,family,genus,species,strain
    Taxonomy = {}
    UniqueTax = []
    sourmashTSV = os.path.join(args.workdir, 'sourmash.csv')
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
    status('Found {:} taxonomic classifications for contigs:\n{:}'.
                format(len(UniqueTax), '\n'.join(UniqueTax)))
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
    status('Dropping {:} contigs from taxonomy screen'.format(len(Tax2Drop)))
    sourTax = os.path.join(args.workdir, 'sourmashed-tax-screen.fasta')
    with open(sourTax, 'w') as outfile:
        with open(os.path.join(args.workdir,assembly_working), 'rU') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                if not record.id in Tax2Drop:
                    SeqIO.write(record, outfile, 'fasta')

    # only do coverage trimming if reads provided
    Contigs2Drop = [] # this will be empty if no reads given to gather by coverage
    if forReads:
        #check if BAM present, if so skip running
        if not os.path.isfile(os.path.join(args.workdir, blobBAM)):
            # index
            bwa_index  = ['bwa','index', os.path.basename(sourTax)]
            status('Building BWA index')
            printCMD(bwa_index)
            subprocess.run(bwa_index, cwd=args.workdir, stderr=DEVNULL)
            #mapped reads to assembly using BWA
            bwa_cmd = ['bwa','mem',
                       '-t', str(args.cpus),
                       os.path.basename(sourTax), # assembly index base
                       forReads]
            if revReads:
                bwa_cmd.append(revReads)

                #run BWA and pipe to samtools sort
                status('Aligning reads to assembly with BWA')
                printCMD(bwa_cmd)
                p1 = subprocess.Popen(bwa_cmd, cwd=args.workdir,
                                      stdout=subprocess.PIPE, stderr=DEVNULL)
                p2 = subprocess.Popen(['samtools', 'sort',
                                       '--threads', str(bamthreads),
                                       '-o', blobBAM, '-'], cwd=args.workdir,
                                      stdout=subprocess.PIPE, stderr=DEVNULL,
                                      stdin=p1.stdout)
                p1.stdout.close()
                p2.communicate()
                subprocess.run(['samtools', 'index', blobBAM], cwd=args.workdir)

        #now calculate coverage from BAM file
        status('Calculating read coverage per contig')
        FastaBed = os.path.join(args.workdir, 'assembly.bed')
        lengths = []
        with open(FastaBed, 'w') as bedout:
            with open(sourTax, 'rU') as SeqIn:
                for record in SeqIO.parse(SeqIn, 'fasta'):
                    bedout.write('{:}\t{:}\t{:}\n'.format(record.id, 0, len(record.seq)))
                    lengths.append(len(record.seq))

        N50 = calcN50(lengths)
        Coverage = {}
        coverageBed = os.path.join(args.workdir, 'coverage.bed')
        cov_cmd = ['samtools', 'bedcov', os.path.basename(FastaBed), blobBAM]
        printCMD(cov_cmd)
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
        # should we make this a variable? 5% was something arbitrary
        min_coverage = float(n50AvgCov * minpct)
        status('Average coverage for N50 contigs is {:}X'.format(int(n50AvgCov)))

        #Start list of contigs to drop
        for k,v in Coverage.items():
            if v[1] <= min_coverage:
                Contigs2Drop.append(k)
        status('Found {:,} contigs with coverage less than {:.2f}X ({:}%)'.
                            format(len(Contigs2Drop), min_coverage, args.mincovpct))

    if args.debug:
        print('Contigs dropped due to coverage: {:}'.format(','.join(Contigs2Drop)))
        print('Contigs dropped due to taxonomy: {:}'.format(','.join(Tax2Drop)))

    DropFinal = Contigs2Drop + Tax2Drop
    DropFinal = set(DropFinal)
    status('Dropping {:,} total contigs based on taxonomy and coverage'.format(len(DropFinal)))
    with open(args.outfile, 'w') as outfile, open(sourTax, 'rU') as seqin:
        for record in SeqIO.parse(seqin, 'fasta'):
            if not record.id in DropFinal:
                SeqIO.write(record, outfile, 'fasta')

    numSeqs, assemblySize = fastastats(args.outfile)
    status('Sourpurged assembly is {:,} contigs and {:,} bp'.
                format(numSeqs, assemblySize))
    if '_' in args.outfile:
        nextOut = args.outfile.split('_')[0]+'.rmdup.fasta'
    elif '.' in args.outfile:
        nextOut = args.outfile.split('.')[0]+'.rmdup.fasta'
    else:
        nextOut = args.outfile+'.rmdup.fasta'

    if checkfile(sourmashTSV):
        baseinput = os.path.basename(args.input)
        if '.' in baseinput:
            baseinput = baseinput.rsplit('.',1)[0]

        shutil.copy(sourmashTSV, baseinput+'.sourmash-taxonomy.csv')

    if not args.debug:
        SafeRemove(args.workdir)

    if not args.pipe:
        status('Your next command might be:\n\tAAFTF rmdup -i {:} -o {:}\n'.format(args.outfile, nextOut))
