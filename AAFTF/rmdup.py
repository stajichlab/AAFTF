import sys
import os
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
import operator


# this runs rountines to identify and remove duplicate
# contigs

from AAFTF.utility import calcN50
from AAFTF.utility import softwrap
from AAFTF.utility import fastastats
from AAFTF.utility import execute
from AAFTF.utility import status
from AAFTF.utility import printCMD
from AAFTF.utility import SafeRemove

def run(parser,args):

    def generateFastas(fasta, pref, query, reference):
        qfile = os.path.join(args.workdir, pref +'query.fasta')
        rfile = os.path.join(args.workdir, pref +'reference.fasta')
        with open(qfile, 'w') as qout:
            with open(rfile, 'w') as rout:
                with open(fasta, 'rU') as infile:
                    for Header,Seq in SimpleFastaParser(infile):
                        if Header in query:
                            qout.write('>{:}\n{:}\n'.format(Header, softwrap(Seq)))
                        elif Header in reference:
                            rout.write('>{:}\n{:}\n'.format(Header, softwrap(Seq)))
        return qfile, rfile
    
    def runMinimap2(query, reference, name):
        FNULL = open(os.devnull, 'w')
        garbage = False #assume this is a good contig
        for line in execute(['minimap2', '-t', str(args.cpus), '-x', 'asm5', '-N5', reference, query], '.'):
            qID, qLen, qStart, qEnd, strand, tID, tLen, tStart, tEnd, matches, alnLen, mapQ = line.split('\t')[:12]
            pident = float(matches) / int(alnLen) * 100
            cov = float(alnLen) / int(qLen) * 100
            if args.debug:
                print('\tquery={:} hit={:} pident={:.2f} coverage={:.2f}'.format(qID, tID, pident, cov))
            
                if pident > args.percent_id and cov > args.percent_cov:
                    print("{:} duplicated: {:.0f}% identity over {:.0f}% of the contig. length={:}".format(name, pident, cov, qLen))
                    garbage = True
                    break
        return garbage #false is good, true is repeat


    #start here -- functions nested so they can inherit the arguments
    custom_workdir = 1
    if not args.workdir:
        custom_workdir = 0
        args.workdir = 'aaftf-rmdup_'+str(os.getpid())
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)
        
    if args.debug:
        status(args)
    status('Looping through assembly shortest --> longest searching for duplicated contigs using minimap2')
    numSeqs, assemblySize = fastastats(args.input)
    fasta_lengths = []
    with open(args.input, 'rU') as infile:
        for Header, Seq in SimpleFastaParser(infile):
            fasta_lengths.append(len(Seq))
    n50 = calcN50(fasta_lengths, num=0.75)
    status('Assembly is {:,} contigs; {:,} bp; and N75 is {:,} bp'.format(numSeqs, assemblySize, n50))
    
    #get list of tuples of sequences sorted by size (shortest --> longest)
    AllSeqs = {}
    with open(args.input, 'rU') as infile:
        for Header, Seq in SimpleFastaParser(infile):
            if not Header in AllSeqs:
                AllSeqs[Header] = len(Seq)
    sortSeqs = sorted(AllSeqs.items(), key=operator.itemgetter(1),reverse=False)
    if args.exhaustive:
        n50 = sortSeqs[-1][1]
    those2check = [x for x in sortSeqs if x[1] < n50]
    status('Will check {:,} contigs for duplication --> those that are < {:,} && > {:,}'.format(len(those2check), n50, args.minlen))
    #loop through sorted list of tuples
    ignore = []
    for i,x in enumerate(sortSeqs):
        sys.stdout.flush()
        if x[1] < args.minlen:
            ignore.append(x[0])
            continue
        if x[1] > n50:
            sys.stdout.flush()
            sys.stdout.write('\n')
            break
        if args.debug:
            status('Working on {:} len={:} remove_tally={:}'.format(x[0], x[1], len(ignore)))
        else:
            text = "\rProgress: {:} of {:}; remove tally={:,}; current={:}; length={:}".format(i, len(those2check), len(ignore), x[0], x[1])
            sys.stdout.write(text)
        #generate input files for minimap2
        theRest = [i[0] for i in sortSeqs[i+1:]]
        pid = str(os.getpid())
        qfile, rfile = generateFastas(args.input, pid, x[0], theRest)
        #run minimap2
        result = runMinimap2(qfile, rfile, x[0])
        if result:
            ignore.append(x[0])
    
    ignore = set(ignore)
    with open(args.out, 'w') as clean_out:
        with open(args.input, 'rU') as infile:
            for Header, Seq in SimpleFastaParser(infile):
                if not Header in ignore:
                    clean_out.write('>{:}\n{:}\n'.format(Header, softwrap(Seq)))
    numSeqs, assemblySize = fastastats(args.out)
    status('Cleaned assembly is {:,} contigs and {:,} bp'.format(numSeqs, assemblySize))
    if '_' in args.out:
        nextOut = args.out.split('_')[0]+'.pilon.fasta'
    elif '.' in args.out:
        nextOut = args.out.split('.')[0]+'.pilon.fasta'
    else:
        nextOut = args.out+'.pilon.fasta'
    status('Your next command might be:\n\tAAFTF pilon -i {:} -l PE_R1.fastq.gz -r PE_R2.fastq.gz -o {:}\n'.format(args.out, nextOut))

    if not args.debug and not custom_workdir:
        SafeRemove(args.workdir)        
