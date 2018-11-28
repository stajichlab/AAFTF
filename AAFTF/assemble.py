# run a set of default genome assembly using
# SPAdes. Additional tools could be supported or
# users may prefer to run their custom assembly and skip this step


import sys, os, subprocess, shutil
from AAFTF.utility import status
from AAFTF.utility import printCMD
from AAFTF.utility import fastastats

def run(parser,args):

    if not args.workdir:
        args.workdir = 'spades_'+str(os.getpid())
        
    spadescmd = ['spades.py','--threads', str(args.cpus), 
                 '--cov-cutoff', 'auto',
                 '--mem', args.memory, '--careful', '-o', args.workdir]

    if args.spades_tmpdir:
        spadescmd.extend(['--tmp-dir',args.spades_tmpdir])
        
    #find reads -- use --left/right or look for cleaned in tmpdir
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        status('Unable to located FASTQ raw reads, provide --left')
        sys.exit(1)
    
    if not revReads:
        spadescmd = spadescmd + ['-s', forReads]
    else:
        spadescmd = spadescmd + ['--pe1-1', forReads, '--pe1-2', revReads]

    if os.path.isdir(args.workdir):
        spadescmd.append('--continue')

    # now run the spades job
    status('Assembling FASTQ data using Spades')
    printCMD(spadescmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(spadescmd)
    else:
        subprocess.run(spadescmd, stdout=DEVNULL, stderr=DEVNULL)
    #pull out assembly
    if args.out:
        finalOut = args.out
    else:
        finalOut = prefix+'.spades.fasta'

    if os.path.isfile(os.path.join(args.workdir, 'scaffolds.fasta')):
        shutil.copyfile(os.path.join(args.workdir,'scaffolds.fasta'), finalOut)
        status('Spades assembly finished: {:}'.format(finalOut))
        numSeqs, assemblySize = fastastats(finalOut)
        status('Assembly is {:,} scaffolds and {:,} bp'.format(numSeqs, assemblySize))
    else:
        status('Spades assembly output missing -- check Spades logfile.')
        
    if not args.pipe:
        status('Your next command might be:\n\tAAFTF vecscreen -i {:} -c {:}\n'.format(finalOut, args.cpus))
    
