"""Run the fcs-gx tool to look for contaminants."""
import os
import shutil
import subprocess
import sys
import uuid

from Bio import SeqIO

from AAFTF.utility import (SafeRemove, calcN50, checkfile, execute, fastastats, printCMD, status)

# logging - we may need to think about whether this has
# separate name for the different runfolder

# flake8: noqa: C901
def run(parser, args):
    """Run NCBI fcs_gx routines to detect and remove contaminant contigs."""
    if not args.workdir:
        args.workdir = f'aaftf-fcsgx_{str(uuid.uuid4())[:8]}'
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    # parse database locations
    if not args.db or not os.path.isfile(f'{args.db}.gxi'):
        status(f"{args.db}.gxi not found, needs to have setup the fcs_gx db - https://github.com/ncbi/fcs/wiki/FCS-GX")
        sys.exit(1)

    numSeqs, assemblySize = fastastats(os.path.join(args.input))
    status('Assembly is {:,} contigs and {:,} bp'.format(numSeqs,assemblySize))
    DEVNULL = open(os.devnull, 'w')

    #now filter for taxonomy with sourmash lca classify
    status('Running fcs_gx to get taxonomy classification for each contig')
    
    #python scripts/run_gx.py --bin-dir dist --gx-db /sw/db/gxdb --tax-id 4842 --fasta

    fcsgx_compute = ['run_gx', '--fasta', args.input, '--tax-id', args.taxid,
                    '--gx-db', args.db, '--out-dir', args.workdir]
    printCMD(fcsgx_compute)
    fcs_log = "fcs_gx.log"
    with open(os.path.join(args.workdir, fcs_log), 'w') as logfile:
        subprocess.run(fcsgx_compute, stderr=logfile)

    fname = os.path.splitext(os.path.basename(args.input))[0]
    # output tsv:
    # seq_id	start_pos	end_pos	seq_len	action	div	agg_cont_cov	top_tax_name

    Seq2Drop = {}

    fcsgxTSV = os.path.join(args.workdir, f'{fname}.{args.taxid}.fcs_gx_report.txt')
    if not os.path.isfile(fcsgxTSV):
        status(f'fcs_gx did not produce file {fcsgxTSV}')
        return
    with open(fcsgxTSV) as fcsgx_out:
        for line in fcsgx_out:
            if not line.startswith("#"):
                line = line.strip()
                cols = line.split("\t")
                Seq2Drop[cols[0]] = 1

    # drop contigs from taxonomy before calculating coverage
    status(f'Dropping {len(Seq2Drop)} contigs from fcs-gx taxonomy screen')
    with open(args.outfile, 'w') as ofh:
        for record in SeqIO.parse(args.input, 'fasta'):
            if record.id not in Seq2Drop:
                SeqIO.write(record, ofh, 'fasta')

    if args.debug:
        print('Contigs dropped due to taxonomy: {:}'.format(','.join(Seq2Drop)))

    numSeqs, assemblySize = fastastats(args.outfile)
    status('fcs-gx assembly is {:,} contigs and {:,} bp'.
                format(numSeqs, assemblySize))
    if '_' in args.outfile:
        nextOut = args.outfile.split('_')[0]+'.rmdup.fasta'
    elif '.' in args.outfile:
        nextOut = args.outfile.split('.')[0]+'.rmdup.fasta'
    else:
        nextOut = f'{args.outfile}.rmdup.fasta'

    if checkfile(fcsgxTSV):
        baseinput = os.path.basename(args.input)
        basedir = os.path.dirname(args.input)
        if '.' in baseinput:
            baseinput = baseinput.rsplit('.',1)[0]

        shutil.copy(fcsgxTSV, os.path.join(basedir,f'{baseinput}.fcs_gx-taxonomy.tsv'))

    if not args.debug:
        SafeRemove(args.workdir)

    if not args.pipe:
        status(f'Your next command might be:\n\tAAFTF rmdup -i {args.outfile} -o {nextOut}\n')
