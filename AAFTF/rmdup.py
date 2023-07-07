"""RMDup module step removes smaller contigs redudnant with larger ones.

This uses minimap to map small contigs against the database of contigs in an
assembly and removes those which are redundant.
"""
import operator
import os
import sys
import uuid

from Bio.SeqIO.FastaIO import SimpleFastaParser

from AAFTF.utility import (SafeRemove, calcN50, execute, fastastats, softwrap,
                           status)


def run(parser, args):
    """Run routines to identify and remove duplicate contigs."""
    def generateFastas(fasta, pref, query, reference):
        qfile = os.path.join(args.workdir, pref + 'query.fasta')
        rfile = os.path.join(args.workdir, pref + 'reference.fasta')
        with open(qfile, 'w') as qout:
            with open(rfile, 'w') as rout:
                with open(fasta) as infile:
                    for Header, Seq in SimpleFastaParser(infile):
                        if Header in query:
                            qout.write(f'>{Header}\n{softwrap(Seq)}\n')
                        elif Header in reference:
                            rout.write(f'>{Header}\n{softwrap(Seq)}\n')
        return qfile, rfile

    def runMinimap2(query, reference, name):
        """Run minimap2 for matching contigs."""
        garbage = False  # assume this is a good contig
        for line in execute(['minimap2', '-t', str(args.cpus), '-x', 'asm5', '-N5', reference, query], '.'):
            qID, qLen, qStart, qEnd, strand, tID, tLen, tStart, tEnd, matches, alnLen, mapQ = line.split('\t')[:12]
            pident = float(matches) / int(alnLen) * 100
            cov = float(alnLen) / int(qLen) * 100
            if args.debug:
                print(f'\tquery={qID} hit={tID} pident={pident:.2f} coverage={cov:.2f}')

                if pident > args.percent_id and cov > args.percent_cov:
                    print(f"{name} duplicated: {pident:.0f}% identity over {cov:.0f}% of the contig. length={qLen}")
                    garbage = True
                    break
        return garbage  # false is good, true is repeat

    # start here -- functions nested so they can inherit the arguments
    custom_workdir = 1
    if not args.workdir:
        custom_workdir = 0
        args.workdir = 'aaftf-rmdup_'+str(uuid.uuid4())[:8]
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    if args.debug:
        status(args)
    status('Looping through assembly shortest --> longest searching for duplicated contigs using minimap2')
    numSeqs, assemblySize = fastastats(args.input)
    fasta_lengths = []
    with open(args.input) as infile:
        # trunk-ignore(ruff/B007)
        for Header, Seq in SimpleFastaParser(infile):
            fasta_lengths.append(len(Seq))
    n50 = calcN50(fasta_lengths, num=0.75)
    status(f'Assembly is {numSeqs:,} contigs; {assemblySize:,} bp; and N75 is {n50:,} bp')

    # get list of tuples of sequences sorted by size (shortest --> longest)
    AllSeqs = {}
    with open(args.input) as infile:
        for Header, Seq in SimpleFastaParser(infile):
            if Header not in AllSeqs:
                AllSeqs[Header] = len(Seq)
    sortSeqs = sorted(AllSeqs.items(), key=operator.itemgetter(1), reverse=False)
    if args.exhaustive:
        n50 = sortSeqs[-1][1]
    those2check = [x for x in sortSeqs if x[1] < n50]
    status(f'Will check {len(those2check):,} contigs for duplication --> those that are < {n50:,} && > {args.minlen:,}')
    # loop through sorted list of tuples
    ignore = []
    for i, x in enumerate(sortSeqs):
        sys.stdout.flush()
        if x[1] < args.minlen:
            ignore.append(x[0])
            continue
        if x[1] > n50:
            sys.stdout.flush()
            sys.stdout.write('\n')
            break
        if args.debug:
            status(f'Working on {x[0]} len={x[1]} remove_tally={len(ignore)}')
        else:
            text = f"\rProgress: {i} of {len(those2check)}; remove tally={len(ignore):,}; current={x[0]}; length={x[1]}     "
            sys.stdout.write(text)
        # generate input files for minimap2
        theRest = [i[0] for i in sortSeqs[i+1:]]
        pid = str(os.getpid())
        qfile, rfile = generateFastas(args.input, pid, x[0], theRest)
        # run minimap2
        result = runMinimap2(qfile, rfile, x[0])
        if result:
            ignore.append(x[0])

    ignore = set(ignore)
    with open(args.out, 'w') as clean_out:
        with open(args.input) as infile:
            for Header, Seq in SimpleFastaParser(infile):
                if Header not in ignore:
                    clean_out.write(f'>{Header}\n{softwrap(Seq)}\n')
    numSeqs, assemblySize = fastastats(args.out)
    status(f'Cleaned assembly is {numSeqs:,} contigs and {assemblySize:,} bp')
    if '_' in args.out:
        nextOut = args.out.split('_')[0]+'.polish.fasta'
    elif '.' in args.out:
        nextOut = args.out.split('.')[0]+'.polish.fasta'
    else:
        nextOut = args.out+'.polish.fasta'

    if not args.pipe:
        status(f'Your next command might be:\n\tAAFTF polish -i {args.out} -l PE_R1.fastq.gz -r PE_R2.fastq.gz -o {nextOut}\n')

    if not args.debug and not custom_workdir:
        SafeRemove(args.workdir)
