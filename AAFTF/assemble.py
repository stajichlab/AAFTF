# run a set of default genome assembly using
# SPAdes. Additional tools could be supported or
# users may prefer to run their custom assembly and skip this step


import sys, os, subprocess, shutil

#logging
import logging
logger = logging.getLogger('AAFTF')

from AAFTF.utility import printCMD
from AAFTF.utility import fastastats

def run(parser,args):

    if args.workdir == 'working_AAFTF' and not args.prefix and not args.left and args.cpus == 1: 
        logger.error(' Please provide either -w,--workdir, -p,--prefix, or --left reads.')
        sys.exit(1)

    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)


    spadescmd = ['spades.py','--threads', str(args.cpus),
                 '-k', '21,33,55,77,99,127','--mem',args.memory,'--careful',
                 '-o', os.path.join(args.workdir, 'spades')]
    if os.path.isdir(os.path.join(args.workdir, 'spades')):
        spadescmd.append('--continue')
        
    #find reads -- use --left/right or look for cleaned in tmpdir
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        for file in os.listdir(args.workdir):
            if '_cleaned' in file and file.endswith('q.gz'):
                if '_1' in file:
                    forReads = os.path.abspath(os.path.join(args.workdir, file))
                if '_2' in file:
                    revReads = os.path.abspath(os.path.join(args.workdir, file))
    if not forReads:
        logger.error('Unable to located FASTQ raw reads, provide correct combination of --prefix, --workdir, or --left')
        sys.exit(1)
    
    if not revReads:
        spadescmd = spadescmd + ['-s', forReads]
    else:
        spadescmd = spadescmd + ['--pe1-1', forReads, '--pe1-2', revReads]

    # now run the spades job
    logger.info('Assembling FASTQ data using Spades')
    logger.info('CMD: {:}'.format(printCMD(spadescmd, 10)))
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(spadescmd)
    else:
        subprocess.run(spadescmd, stdout=DEVNULL, stderr=DEVNULL)
    #pull out assembly
    if args.out:
        finalOut = args.out
    else:
        if args.prefix:
            finalOut = args.prefix+'.spades.fasta'
        else:
            finalOut = os.path.basename(forReads).split('.')[0]+'.spades.fasta'
    if os.path.isfile(os.path.join(args.workdir, 'spades', 'scaffolds.fasta')):
        shutil.copyfile(os.path.join(args.workdir, 'spades', 'scaffolds.fasta'), finalOut)
        logger.info(' Spades assembly finished: {:}'.format(finalOut))
        numSeqs, assemblySize = fastastats(finalOut)
        logger.info('Assembly is {:,} contigs and {:,} bp'.format(numSeqs, assemblySize))
    else:
        logger.error(' Spades assembly output missing -- check Spades logfile.')
    logger.info('Your next command might be:\n\tAAFTF vecscreen -i {:} -c {:}\n'.format(finalOut, args.cpus))
    
