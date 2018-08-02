import sys, csv, re, operator,os, subprocess

# biopython needed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# for processing identify vector sequence
default_percent_id_cutoff = 98

#logging
import logging
logger = logging.getLogger('AAFTF')

# this code is based on this post on stackexchange
# https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def merge_intervals(intervals):
    """
    A simple algorithm can be used:
    1. Sort the intervals in increasing order
    2. Push the first interval on the stack
    3. Iterate through intervals and for each one compare current interval
       with the top of the stack and:
       A. If current interval does not overlap, push on to stack
       B. If current interval does overlap, merge both intervals in to one
          and push on to stack
    4. At the end return stack
    """
    merged = []
    sorted_by_lower_bound = sorted(intervals, key=operator.itemgetter(0))

    if not sorted_by_lower_bound:  # no intervals to merge
        return

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            #added this condition branch
            if higher[0] - lower[1] == 1:
                merged[-1] = (lower[0], higher[1])  # replace by merged interval
            #end addition. Also changed if below to elif
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            elif higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged

def parse_clean_blastn(fastafile,blastn,pid_cutoff):

    fastafile = sys.argv[1]
    blastn = sys.argv[2]
    prefix = os.path.splitext(fastafile)[0]
    cleaned = prefix + ".clean.fsa"
    logging = prefix + ".parse.log"

    excludes = {}

    with open(blastn,"r") as vectab:
        rdr = csv.reader(vectab,delimiter="\t")
        for row in rdr:
            idin = row[0]
            pid = row[2]
            if float(pid) < pid_cutoff:
                continue
                # skip this line since PID is lower than needed

            loc = [int(row[6]), int(row[7])]
            if loc[0] > loc[1]:
                loc = [loc[1],loc[0]]

            if idin not in excludes:
                excludes[idin] = []

            excludes[idin].append(loc)

    with open(cleaned, "w") as output_handle, open(logging,"w") as log:
        for record in SeqIO.parse(fastafile, "fasta"):
            record.description = ""
            if record.id in excludes:
                trimloc = excludes[record.id]

            if len(trimloc) > 1:
                trimloc = sorted(merge_intervals(trimloc),
                                 reverse=True,
                                 key=lambda locitem: locitem[0])

            seqlen = len(record)
            for loc in trimloc:
                left = int( loc[0] ) - 1
                right = int( loc[1] )
                newrecord = Seq("",generic_dna)
                log.write("trimming %d to %d in %s len=%d"
                      % (left,right,record.id,len(record)))
                if left == 0:
                    newrecord = record[right-1:]
                elif right == len(record):
                    newrecord = record[:left]
                else:
                    # internal slicing
                    log.write("-->internal slicing :%d .. %d:" % (left,right-1))
                    log.write('  left string is %s' % record[0:left])
                    log.write('  right string is %s' % record[right-1:])
                    newrecord = record[0:left] + record[right-1:]
                    newrecord.id = record.id
                    log.write(" -- new len for %s is %d: %s" % (newrecord.id, len(newrecord),newrecord))
                record = newrecord

            if(len(record) >= 200):
                SeqIO.write(record, output_handle, "fasta")

def make_blastdb(type,file,name):
    indexfile = file
    if type eq 'nucl':
        indexfile += ".pin"
    else:
        indexfile += ".pin"
        
    if ( not os.path.exists(indexfile) or
         os.path.getctime(indexfile) < os.path.getctime(file)):
        subprocess.call(["makeblastdb",'-dbtype',type,
                         '-in',file,'-out',name])

        
def run(parser,args):

    print(args.tmpdir)
    
    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)

    if args.percent_id:
        percentid_cutoff = args.percent_id
    else:
        percentid_cutoff = default_percent_id_cutoff
        
    # Common Euk/Prot contaminats for blastable DB later on
    makeblastdblist = []
    for db,url in DB_Links:
        dbname = os.path.basename(str(url))
        file = os.path.join(args.tmpdir,dbname)
        if file.endswith(".gz"):
            nogz = os.path.splitext(file)[0]
            urllib.request.urlretrieve(url,file)
            subprocess.call(['gunzip',file])
            make_blastdb('nucl',nogz,os.path.join(args.tmpdir,db))
        else:
            if not os.path.exists(file):
                urllib.request.urlretrieve(url,file)
                make_blastdb('nucl',file,os.path.join(args.tmpdir,db))

    prefix =os.basename(args.infile)
    prefix = os.path.splitext(prefix)[0]
    infile = args.infile
    round = 0
    while 1:
        report = os.path.join(args.tmpdir,"%s.r%d.vecscreen.tab"%(prefix,round))
        subprocess.call(['blastn','-task','blastn',
                        '-reward','1','-penalty','-5','-gapopen','3',
                        '-gapextend', '3', '-dust','yes','-soft_masking','true',
                        '-evalue', '700','-searchsp','1750000000000',
                        '-db', os.path.join(args.tmpdir,'UniVec'),
                        '-outfmt', '6', '-num_threads',args.cpus,
                        '-query', infile,
                        '-out', report])
        # this needs to know/return the new fasta file?
        count = parse_clean_blastn(infile,report,percentid_cutoff)
        if count == 0: # if there are no vector matches < than the pid cutoff
            break
        else:
            round += 1
        
# probably run this in python rather than shell out again and depend on gawk but it is fast either way I am guessing
#if [ -f $PREF.r$ROUND.vecscreen.tab ]; then
#	    COUNT=$(awk -v cutoff=$PID_CUTOFF 'BEGIN { sum = 0 } { if ( $3 >= cutoff ) sum +=1 } END { print sum }' $PREF.r$ROUND.vecscreen.tab)
#	fi

#if [[ $COUNT == 0 ]]; then
#	    RUN=0
#	else
#	    NEXTROUND=ROUND+1
#	    parse_clean_blastn(FASTA,BLAST_TAB)
#	    ln -s $PREF.r$ROUND.clean.fsa $PREF.r$NEXTROUND.fna
#	    ROUND=NEXTROUND
#	    FASTA=$PREF.r$ROUND.fna
#	fi
