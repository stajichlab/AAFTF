# this runs rountines to identify contaminant contigs
# which are presumably sequences that were not screened out
# in the filter step.
# vector library UniVec and known or user specified contaminanting
# sequences (by GenBank accession number) can be provided for
# additional cleanup

# default lirbaries screen are located in resources.py
# and include common Euk, Prok, and MITO contaminants

import sys, csv, re, operator, os, gzip
import shutil

from subprocess import call, Popen, PIPE, STDOUT

import urllib.request
from AAFTF.resources import SeqDBs
from AAFTF.resources import DB_Links
from AAFTF.utility import status
from AAFTF.utility import printCMD
from AAFTF.utility import softwrap
from AAFTF.utility import countfasta
from AAFTF.utility import SafeRemove

# biopython needed
from Bio import SeqIO


BlastPercent_ID_ContamMatch = "90.0"
BlastPercent_ID_MitoMatch   = "98.6"
DEVNULL = open(os.devnull, 'w')

#VecScreen matches
'''
Strong Match to Vector
(Expect 1 random match in 1,000,000 queries of length 350 kb.)
Terminal match with Score ≥ 24.
Internal match with Score ≥ 30.
Moderate Match to Vector
(Expect 1 random match in 1,000 queries of length 350 kb.)
Terminal match with Score 19 to 23.
Internal match with Score 25 to 29.
Weak Match to Vector
(Expect 1 random match in 40 queries of length 350 kb.)
Terminal match with Score 16 to 18.
Internal match with Score 23 to 24.
Segment of Suspect Origin
Any segment of fewer than 50 bases between two vector matches or between a match and an end.
'''

def group(lst, n):
    for i in range(0, len(lst), n):
        val = lst[i:i+n]
        if len(val) == n:
            yield tuple(val)

def parse_clean_blastn(fastafile, prefix, blastn, stringent):
    '''
    Blast header rows:
    qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue score qlen
    '''

    cleaned = prefix + ".clean.fsa"
    logging = prefix + ".parse.log"

    excludes = {}
    VecHits = {}
    found_vector_seq = 0
    with open(blastn,"r") as vectab:
        rdr = csv.reader(vectab,delimiter="\t")
        for row in rdr:
            qaccver,saccver,pid,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,score,qlen = row
            if qaccver in contigs_to_remove:
                continue
            #vecscreen https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Moderate
            #says to use score here (I'm interpret as score not bitscore)
            #need to determine if match is terminal or if internal
            loc = [int(qstart), int(qend)]
            if loc[0] > loc[1]:
                loc = [loc[1],loc[0]]
            #check for location
            terminal = False
            position = None
            if loc[0] <= 25:
                terminal = True
                position = '5'
            if (int(qlen) - loc[1]) <= 25:
                terminal = True
                position = '3'
            Match = 0 # weak=0, moderate=1, strong=2
            score = int(score)
            if terminal:
                if score >= 19:
                    Match = 1
                if score >= 24:
                    Match = 2
            else:
                if score >= 25:
                    Match = 1
                if score >= 30:
                    Match = 2
            if Match == 0:
                continue
            if stringent == 'high':
                if Match > 0:
                    found_vector_seq += 1
                    if not qaccver in VecHits:
                        VecHits[qaccver] = [(saccver, int(qlen), loc, int(score), terminal, position)]
                    else:
                        VecHits[qaccver].append((saccver, int(qlen), loc, int(score), terminal, position))
            else:
                if Match > 1:
                    found_vector_seq += 1
                    if not qaccver in VecHits:
                        VecHits[qaccver] = [(saccver, int(qlen), loc, int(score), terminal, position)]
                    else:
                        VecHits[qaccver].append((saccver, int(qlen), loc, int(score), terminal, position))          
    
    trimTerminal = 0
    splitContig = 0
    with open(cleaned, "w") as output_handle, open(logging, "w") as log:
        for record in SeqIO.parse(fastafile, "fasta"):
            FiveEnd = 0
            ThreeEnd = len(record.seq)
            internals = []
            slicer = []
            sInt = []
            Seq = str(record.seq)
            if not record.id in VecHits:
                if len(record.seq) >= 200:
                    output_handle.write('>{:}\n{:}\n'.format(record.id, softwrap(Seq)))
            else:
                #VecHits contains list of tuples of information, if terminal, then just truncate
                #off the closest side. Also, need to check if multiple intervals are within 50
                #bp of each other, that whole interval is removed.
                #should be able to accomplish above with the several rounds that it runs with, 
                #so split on internal and trim terminal. done.
                for hit in VecHits[record.id]:
                    ID,length,loc,score,terminal,pos = hit
                    if terminal and pos == '5':
                        if loc[1] > FiveEnd:
                            FiveEnd = loc[1]
                    elif terminal and pos == '3':
                        if loc[0] < ThreeEnd:
                            ThreeEnd = loc[0]
                    else: #internal hits to add to list
                        if not loc in internals:
                            internals.append(loc)
                #now sort intervals
                sInt = sorted(internals, key=lambda x: int(x[0]))
                #now construct slicing list
                if len(sInt) < 1:
                    slicer = [FiveEnd,ThreeEnd]
                else:
                    slicer = [FiveEnd]
                    for x in sInt:
                        slicer = slicer + x
                    slicer.append(ThreeEnd)
                paired_slicer = list(group(slicer, 2))
                if len(paired_slicer) < 2:
                    status('Terminal trimming {:} to {:}'.format(record.id, paired_slicer))
                    newSeq = Seq[paired_slicer[0][0]:paired_slicer[0][1]]
                    if len(newSeq) >= 200:
                        output_handle.write('>{:}\n{:}\n'.format(record.id, softwrap(newSeq)))
                else:
                    status('Spliting contig {:} into {:}'.format(record.id, paired_slicer))
                    for num,y in enumerate(paired_slicer):
                        newSeq = Seq[y[0]:y[1]]
                        if len(newSeq) >= 200:
                            output_handle.write('>split{:}_{:}\n{:}\n'.format(num+1, record.id, softwrap(newSeq)))

    return (found_vector_seq, cleaned)

def make_blastdb(type,file,name):
    indexfile = name
    if type == 'nucl':
        indexfile += ".nin"
    else:
        indexfile += ".pin"
    
    if not os.path.exists(indexfile) or os.path.getctime(indexfile) < os.path.getctime(file):
        cmd = ['makeblastdb','-dbtype',type,'-in',file,'-out',name]
        printCMD(cmd)
        call(cmd, stdout=DEVNULL, stderr=DEVNULL)

        
def run(parser,args):
    if not args.workdir:
        args.workdir = 'aaftf-vecscreen_'+str(os.getpid())
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)
    
    #parse database locations
    DB = None
    if not args.AAFTF_DB:
        try:
            DB = os.environ["AAFTF_DB"]
        except KeyError:
            if args.AAFTF_DB:
                DB = args.AAFTF_DB
            else:
                pass
    else:
        DB = args.AAFTF_DB
    
    if args.percent_id:
        percentid_cutoff = args.percent_id

    infile = args.infile
    outfile = os.path.basename(args.outfile)
    outdir = os.path.dirname(args.outfile)
    if '.f' in outfile:
        prefix = outfile.rsplit('.f', 1)[0]
        print("prefix is ",prefix)
    else:
        prefix = str(os.getpid())
    if not outfile:
        outfile = "%s.vecscreen.fasta" % prefix

    outfile_vec = os.path.join(args.workdir,
                               "%s.tmp_vecscreen.fasta" % (prefix))

    # Common Euk/Prot contaminats for blastable DB later on
    status('Building BLAST databases for contamination screen.')
    makeblastdblist = []
    for d in DB_Links:
        if d == 'sourmash':
            continue
        url = DB_Links[d]
        dbname = os.path.basename(str(url))
        #logger.debug("testing for url=%s dbname=%s"%(url,dbname))
        if DB:
            file = os.path.join(DB, dbname)
        else:
            file = os.path.join(args.workdir,dbname)
        if file.endswith(".gz"):
            nogz = os.path.splitext(file)[0]
            if not os.path.exists(nogz):
                if not os.path.exists(file):
                    urllib.request.urlretrieve(url,file)
                
                with gzip.open(file, 'rb') as ingz, open(nogz,'wb') as outfa:
                    shutil.copyfileobj(ingz,outfa)
#                call(['gunzip', '-k', file])
                make_blastdb('nucl', nogz, os.path.join(args.workdir,d))
            else:
                make_blastdb('nucl', nogz, os.path.join(args.workdir,d))
        else:
            if not os.path.exists(file):
                urllib.request.urlretrieve(url,file)
            make_blastdb('nucl',file,os.path.join(args.workdir,d))
    
    global contigs_to_remove
    contigs_to_remove = {}
    regions_to_trim = {}
    
    #qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
    for contam in ["CONTAM_EUKS","CONTAM_PROKS" ]:                       
        status("%s Contamination Screen" % (contam))
        blastreport = os.path.join(args.workdir,
                                   "%s.%s.blastn" % (contam, prefix))
        blastnargs = ['blastn',
                      '-query', infile,
                      '-db', os.path.join(args.workdir,contam),
                      '-num_threads', str(args.cpus),
                      '-dust', 'yes', '-soft_masking', 'true',
                      '-perc_identity',BlastPercent_ID_ContamMatch,
                      '-lcase_masking', '-outfmt', '6', '-out',blastreport]
        printCMD(blastnargs)
        call(blastnargs)
        hits = 0
        with open(blastreport) as report:
            colparser = csv.reader(report, delimiter="\t")
            for row in colparser:
                if( ( float(row[2]) >= 98.0 and
                      int(row[3]) >= 50)  or
                    ( float(row[2]) >= 94.0 and
                      int(row[3]) >= 100) or
                    ( float(row[2]) >= 90.0 and
                      int(row[3]) >= 200) ):
                    if not row[0] in regions_to_trim:
                        if int(row[6]) < int(row[7]):
                            start = int(row[6])
                            end = int(row[7])
                        else:
                            start = int(row[7])
                            end = int(row[6])
                        regions_to_trim[row[0]] = [(start, end, contam, row[1], float(row[2]))]
                    else:
                        regions_to_trim[row[0]].append((start, end, contam, row[1], float(row[2])))
        status('{:} screening finished'.format(contam))

    eukCleaned = os.path.join(args.workdir, "%s.euk-prot_cleaned.fasta" % (prefix))
    if len(regions_to_trim) > 0:
        with open(eukCleaned, 'w') as cleanout:
            with open(infile, 'rU') as fastain:
                for record in SeqIO.parse(fastain, 'fasta'):
                    if not record.id in regions_to_trim:
                        cleanout.write('>{:}\n{:}\n'.format(record.id, softwrap(str(record.seq))))
                    else:
                        Seq = str(record.seq)
                        regions = regions_to_trim[record.id]
                        status('Splitting {:} due to contamination: {:}'.format(record.id, regions))
                        lastpos = 0
                        newSeq = ''
                        for i,x in enumerate(regions):
                            newSeq = Seq[lastpos:x[0]]
                            lastpos = x[1]
                            cleanout.write('>split{:}_{:}\n{:}\n'.format(i, record.id, softwrap(newSeq)))
                            if i == len(regions)-1:
                                newSeq = Seq[x[1]:]
                                cleanout.write('>split{:}_{:}\n{:}\n'.format(i+1, record.id, softwrap(newSeq)))
    else:
        eukCleaned = infile
            
    # MITO screen
    status('Mitochondria Contamination Screen')
    mitoHits = []
    blastreport = os.path.join(args.workdir,
                               "%s.%s.blastn" % ('MITO',prefix))
    blastnargs = ['blastn',
                  '-query', eukCleaned,
                  '-db', os.path.join(args.workdir,'MITO'),
                  '-num_threads', str(args.cpus),
                  '-dust', 'yes', '-soft_masking', 'true',
                  '-perc_identity',BlastPercent_ID_MitoMatch,
                  '-lcase_masking', '-outfmt','6',
                  '-out', blastreport]
    printCMD(blastnargs)
    call(blastnargs)
    with open(blastreport) as report:
        colparser = csv.reader(report, delimiter="\t")
        for row in colparser:
            if int(row[3]) >= 120:
                contigs_to_remove[row[0]] = ('MitoScreen', row[1], float(row[2]))
                mitoHits.append(row[0])
    status('Mito screening finished.')

    #vecscreen starts here
    status('Starting VecScreen, will remove terminal matches and split internal matches')
    rnd = 0
    count = 1
    while (count > 0):
        filepref = "%s.r%d" % (prefix,rnd)
        report = os.path.join(args.workdir,"%s.vecscreen.tab"%(filepref))
        if not os.path.exists(report):
            cmd = ['blastn','-task','blastn',
                  '-reward','1','-penalty','-5','-gapopen','3',
                  '-gapextend', '3', '-dust','yes','-soft_masking','true',
                  '-evalue', '700','-searchsp','1750000000000',
                  '-db', os.path.join(args.workdir,'UniVec'),
                  '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore score qlen', 
                  '-num_threads',str(args.cpus),
                  '-query', eukCleaned, '-out', report]
            #logger.info('CMD: {:}'.format(printCMD(cmd,7)))
            call(cmd)
        # this needs to know/return the new fasta file?
        status("Parsing VecScreen round {:}: {:} for {:}".format(rnd+1, filepref,report))
        (count, cleanfile) = parse_clean_blastn(eukCleaned, os.path.join(args.workdir,filepref),report, args.stringency)
        status("count is %d cleanfile is %s"%(count, cleanfile))
        if count == 0: # if there are no vector matches < than the pid cutoff
            status("copying %s to %s"%(eukCleaned, outfile_vec))
            shutil.copy(eukCleaned, outfile_vec)
        else:
            rnd += 1
            eukCleaned = cleanfile

    status("{:,} contigs will be removed:".format(len(contigs_to_remove)))
    for k,v in sorted(contigs_to_remove.items()):
        print('\t{:} --> dbhit={:}; hit={:}; pident={:}'.format(k,v[0], v[1], v[2]))

    # this could instead use the outfile and strip .fasta/fsa/fna and add mito on it I suppose, but assumes 
    # a bit about the naming structure

    mitochondria = os.path.join(outdir,prefix+'.mitochondria.fasta')
    with open(args.outfile, "w") as output_handle, open(mitochondria, 'w') as mito_handle:
        for record in SeqIO.parse(outfile_vec, "fasta"):
            if not record.id in contigs_to_remove:
                SeqIO.write(record, output_handle, "fasta")
            elif record.id in mitoHits:
                SeqIO.write(record, mito_handle, "fasta")
    status('Writing {:,} cleaned contigs to: {:}'.format(countfasta(args.outfile), args.outfile))
    status('Writing {:,} mitochondrial contigs to: {:}'.format(countfasta(mitochondria), mitochondria))
    if '_' in args.outfile:
        nextOut = args.outfile.split('_')[0]+'.sourpurge.fasta'
    elif '.' in args.outfile:
        nextOut = args.outfile.split('.')[0]+'.sourpurge.fasta'
    else:
        nextOut = args.outfile+'.sourpurge.fasta'

    if not args.pipe:
        status('Your next command might be:\n\tAAFTF sourpurge -i {:} -o {:} -c {:} --phylum Ascomycota\n'.format(
            args.outfile, nextOut, args.cpus))
    
    if not args.debug:
        SafeRemove(args.workdir)
