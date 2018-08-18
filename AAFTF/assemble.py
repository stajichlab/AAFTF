# run a set of default genome assembly using
# SPAdes. Additional tools could be supported or
# users may prefer to run their custom assembly and skip this step


import sys, os, subprocess, shutil

#logging
import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):
    
    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)


    spadescmd = ['spades.py','--threads', str(args.cpus),
                 '-k', '21,33,55,77,99,127','--mem',args.memory,'--careful',
                 '-o', os.path.join(args.tmpdir, 'spades')]
    if os.path.isdir(os.path.join(args.tmpdir, 'spades')):
        spadescmd.append('--continue')
        
    #find reads -- use --left/right or look for cleaned in tmpdir
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        for file in os.listdir(args.tmpdir):
            if '_cleaned' in file and file.endswith('q.gz'):
                if '_1' in file:
                    forReads = os.path.abspath(os.path.join(args.tmpdir, file))
                if '_2' in file:
                    revReads = os.path.abspath(os.path.join(args.tmpdir, file))
    if not forReads:
        print('Unable to located FASTQ raw reads')
        sys.exit(1)
    
    if not revReads:
        spadescmd = spadescmd + ['-s', forReads]
    else:
        spadescmd = spadescmd + ['--pe1-1', forReads, '--pe1-2', revReads]

    # now run the spades job
    logger.info('CMD: {:}'.format(' '.join(spadescmd)))
    subprocess.run(spadescmd)
    
    #pull out assembly
    if args.out:
        finalOut = args.out
    else:
        if args.prefix:
            finalOut = args.prefix+'.spades.fasta'
        else:
            finalOut = os.basename(forReads).split('.')[0]+'.spades.fasta'
    if os.path.isfile(os.path.join(args.tmpdir, 'spades', 'scaffolds.fasta')):
        shutil.copyfile(os.path.join(args.tmpdir, 'spades', 'scaffolds.fasta'), finalOut)
        logger.info(' Spades assembly finished: {:}'.format(finalOut))
    else:
        logger.error(' Spades assembly output missing -- check logfiles.')
    
