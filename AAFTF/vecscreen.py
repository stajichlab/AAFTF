# this runs rountines to identify contaminant contigs
# which are presumably sequences that were not screened out
# in the filter step.
# vector library UniVec and known or user specified contaminanting
# sequences (by GenBank accession number) can be provided for
# additional cleanup

# default lirbaries screen are located in resources.py
# and include common Euk, Prok, and MITO contaminants

import sys, csv, re, operator,os
import shutil

from subprocess import call, Popen, PIPE, STDOUT

import urllib.request
from AAFTF.resources import SeqDBs
from AAFTF.resources import DB_Links

# biopython needed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# for processing identify vector sequence
default_percent_id_cutoff = 98

#logging
import logging
logger = logging.getLogger('AAFTF')

BlastPercent_ID_ContamMatch = "90.0"
BlastPercent_ID_MitoMatch   = "98.6"
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

def parse_clean_blastn(fastafile,prefix,blastn,pid_cutoff):

    cleaned = prefix + ".clean.fsa"
    logging = prefix + ".parse.log"

    excludes = {}
    found_vector_seq = 0
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
            found_vector_seq += 1

    if found_vector_seq == 0:
        return (0,fastafile)

    with open(cleaned, "w") as output_handle, open(logging,"w") as log:
        for record in SeqIO.parse(fastafile, "fasta"):
            record.description = ""
            trimloc = []
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

    return (found_vector_seq,cleaned)

def make_blastdb(type,file,name):
    indexfile = name
    if type == 'nucl':
        indexfile += ".nin"
    else:
        indexfile += ".pin"
    
    if ( not os.path.exists(indexfile) or
         os.path.getctime(indexfile) < os.path.getctime(file)):
        call(["makeblastdb",'-dbtype',type,
              '-in',file,'-out',name])

        
def run(parser,args):

    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    
    if args.percent_id:
        percentid_cutoff = args.percent_id
    else:
        percentid_cutoff = default_percent_id_cutoff
        
    prefix = os.path.basename(args.infile)
    prefix = os.path.splitext(prefix)[0]
    infile = args.infile
    outfile = args.outfile

    if not outfile:
        outfile = "%s.vecscreen.fasta" % prefix

    outfile_vec = os.path.join(args.tmpdir,
                               "%s.tmp_vecscreen.fasta" % (prefix))

    # Common Euk/Prot contaminats for blastable DB later on
    makeblastdblist = []
    for db in DB_Links:
        url = DB_Links[db]
        dbname = os.path.basename(str(url))
        file = os.path.join(args.tmpdir,dbname)
        if file.endswith(".gz"):
            nogz = os.path.splitext(file)[0]
            if not os.path.exists(nogz):
                if not os.path.exists(file):
                    urllib.request.urlretrieve(url,file)
                call(['gunzip',file])
            make_blastdb('nucl',nogz,os.path.join(args.tmpdir,db))
        else:
            if not os.path.exists(file):
                urllib.request.urlretrieve(url,file)
            make_blastdb('nucl',file,os.path.join(args.tmpdir,db))

    rnd = 0
    count = 1
    while (count > 0):
        filepref = "%s.r%d" % (prefix,rnd)
        report = os.path.join(args.tmpdir,"%s.vecscreen.tab"%(filepref))
        if not os.path.exists(report):
            call(['blastn','-task','blastn',
                  '-reward','1','-penalty','-5','-gapopen','3',
                  '-gapextend', '3', '-dust','yes','-soft_masking','true',
                  '-evalue', '700','-searchsp','1750000000000',
                  '-db', os.path.join(args.tmpdir,'UniVec'),
                  '-outfmt', '6', '-num_threads',str(args.cpus),
                  '-query', infile, '-out', report])
        # this needs to know/return the new fasta file?
        logger.info("parsing %s for %s infile=%s"%(filepref,report,infile))
        (count,cleanfile) = parse_clean_blastn(infile,
                                               os.path.join(args.tmpdir,
                                                            filepref),
                                               report,percentid_cutoff)       
        logger.info("count is %d cleanfile is %s"%(count,cleanfile))
        if count == 0: # if there are no vector matches < than the pid cutoff
            logger.info("copying %s to %s"%(infile,outfile_vec))
            shutil.copy(infile,outfile_vec)
        else:
            rnd += 1
            infile = cleanfile

    # loop is finished for vecscreen now screen for common contaminants

    contigs_to_remove = {}
    for contam in ["CONTAM_EUKS",
                   "CONTAM_PROKS" ]:                       
        logger.info("%s Contamination Screen" % (contam))
        blastreport = os.path.join(args.tmpdir,
                                   "%s.%s.blastn" % (contam, prefix))
        blastnargs = ['blastn',
                      '-query', outfile_vec,
                      '-db', os.path.join(args.tmpdir,contam),
                      '-num_threads', str(args.cpus),
                      '-dust', 'yes', '-soft_masking', 'true',
                      '-perc_identity',BlastPercent_ID_ContamMatch,
                      '-lcase_masking', '-outfmt',
                      '6 "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                      '-out',blastreport ]
        call(blastnargs)
        
        with open(blastreport) as report:
            colparser = csv.reader(report, delimiter="\t")
            for row in colparser:
                if( ( float(row[2]) >= 98.0 and
                      int(row[3]) >= 50)  or
                    ( float(row[2]) >= 94.0 and
                      int(row[3]) >= 100) or
                    ( float(row[2]) >= 90.0 and
                      row[3] >= 200) ):
                    logger.info("Removing contig %s because of match to %s"%(row[0],row[1]))
                    contigs_to_remove[row[0]] = 1
                    
            #done with EUK and PROK screen
            
    # MITO screen
    blastreport = os.path.join(args.tmpdir,
                               "%s.%s.blastn" % ('MITO',prefix))
    blastnargs = ['blastn',
                  '-query', outfile_vec,
                  '-db', os.path.join(args.tmpdir,'MITO'),
                  '-num_threads', str(args.cpus),
                  '-dust', 'yes', '-soft_masking', 'true',
                  '-perc_identity',BlastPercent_ID_MitoMatch,
                  '-lcase_masking', '-outfmt','6',
                  '-out', blastreport]
    call(blastnargs)
    with open(blastreport) as report:
        colparser = csv.reader(report, delimiter="\t")
        for row in colparser:
            if int(row[3]) >= 120:
                logger.info("Removing contig %s because of match to %s"%(row[0],row[1]))
                contigs_to_remove[row[0]] = 1
        
    print("These contigs will be skipped:",contigs_to_remove)

    with open(outfile, "w") as output_handle:
        for record in SeqIO.parse(outfile_vec, "fasta"):
            if record.id not in contigs_to_remove:
                SeqIO.write(record, output_handle, "fasta")
