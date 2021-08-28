# run a set of default genome assembly using
# SPAdes. Additional tools could be supported or
# users may prefer to run their custom assembly and skip this step


import sys, os, subprocess, shutil, uuid
from AAFTF.utility import status
from AAFTF.utility import printCMD
from AAFTF.utility import fastastats

# run spades
def run_spades(parser,args):

    if not args.workdir:
        args.workdir = 'spades_'+str(uuid.uuid4())[:8]

    runcmd = ['spades.py','--threads', str(args.cpus), '--mem', args.memory,
              '-o',args.workdir]

    if args.assembler_args:
        runcmd.extend(args.assembler_args)

    if '--meta' not in runcmd:
        runcmd.extend(['--cov-cutoff','auto',  '--careful'])

    if args.tmpdir:
        runcmd.extend(['--tmp-dir',args.tmpdir])

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
        runcmd = runcmd + ['-s', forReads]
    else:
        runcmd = runcmd + ['--pe1-1', forReads, '--pe1-2', revReads]

    # this basically overrides everything above and only runs --restart-from option
    if os.path.isdir(args.workdir):
        runcmd = [ 'spades.py', '-o', args.workdir,
                   '--threads', str(args.cpus),
                   '--mem', args.memory,
                   '--restart-from last' ]

    # now run the spades job
    status('Assembling FASTQ data using Spades')
    printCMD(runcmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(runcmd)
    else:
        subprocess.run(runcmd, stdout=DEVNULL, stderr=DEVNULL)
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



# run dipspades
def run_dipspades(parser,args):

    if not args.workdir:
        args.workdir = 'dipspades_'+str(os.getpid())

    runcmd = ['dipspades.py','--threads', str(args.cpus), '--cov-cutoff',
              'auto', '--mem', args.memory, '-o',
              args.workdir]

    if args.assembler_args:
        runcmd.extend(args.assembler_args)

    if args.haplocontigs:
        runcmd.extend(['--hap',args.haplocontigs])

    if args.tmpdir:
        runcmd.extend(['--tmp-dir',args.tmpdir])

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
        runcmd = runcmd + ['-s', forReads]
    else:
        runcmd = runcmd + ['--pe1-1', forReads, '--pe1-2', revReads]

        # this basically overrides everything above and only runs --restart-from option
    if os.path.isdir(args.workdir):
        runcmd = [ 'dipspades.py', '-o', args.workdir,
                   '--continue']

    # now run the spades job
    status('Assembling FASTQ data using Spades')

    printCMD(runcmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(runcmd)
    else:
        subprocess.run(runcmd, stdout=DEVNULL, stderr=DEVNULL)
    #pull out assembly

    if args.out:
        finalOut = args.out
    else:
        finalOut = prefix+'.dipspades.fasta'
    dipspadesoutdir = os.path.join(args.workdir,'dipspades')
    if os.path.isfile(os.path.join(args.workdir, 'consensus_contigs.fasta')):
        shutil.copyfile(os.path.join(args.workdir,'consensus_contigs.fasta'), finalOut)
        shutil.copyfile(os.path.join(args.workdir,'dipspades','paired_consensus_contigs.fasta'), prefix+".dipspades_consensus_paired.fasta")
        shutil.copyfile(os.path.join(args.workdir,'dipspades','paired_consensus_contigs.fasta'), prefix+".dipspades_consensus_unpaired.fasta")
        status('Dipspades assembly finished: {:}'.format(finalOut))
        status('Dipspades assembly copied over: {:}'.format(prefix+".dipspades_consensus_unpaired.fasta"),prefix+".dipspades_consensus_paired.fasta")
        numSeqs, assemblySize = fastastats(finalOut)
        status('Assembly is {:,} scaffolds and {:,} bp'.format(numSeqs, assemblySize))
    else:
        status('Spades assembly output missing -- check Dipspades logfile in {:}.'.format(os.path.join(args.workdir,'dipspades','dipspades.log')))

    if not args.pipe:
        status('Your next command might be:\n\tAAFTF vecscreen -i {:} -c {:}\n'.format(finalOut, args.cpus))

def run_megahit(parser,args):

    if not args.workdir:
        args.workdir = 'megahit_'+str(os.getpid())

    runcmd = ['megahit','-t', str(args.cpus),
              '-o', args.workdir]

    if args.assembler_args:
        runcmd.extend(args.assembler_args)

    if args.memory:
        runcmd.extend(['--memory',args.memory])

    if args.tmpdir:
        runcmd.extend(['--tmp-dir',args.tmpdir])

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
        runcmd = runcmd + ['-r', forReads]
    else:
        runcmd = runcmd + ['-1', forReads, '-2', revReads]

    if os.path.isdir(args.workdir):
        status("Cannot re-run with existing folder {}".format(args.workdir))

    # now run the spades job
    status('Assembling FASTQ data using megahit')
    printCMD(runcmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(runcmd)
    else:
        subprocess.run(runcmd, stdout=DEVNULL, stderr=DEVNULL)
    #pull out assembly
    if args.out:
        finalOut = args.out
    else:
        finalOut = prefix+'.megahit.fasta'

    if os.path.isfile(os.path.join(args.workdir, 'final.contigs.fa')):
        shutil.copyfile(os.path.join(args.workdir,'final.contigs.fa'), finalOut)
        status('Megahit assembly finished: {:}'.format(finalOut))
        numSeqs, assemblySize = fastastats(finalOut)
        status('Assembly is {:,} scaffolds and {:,} bp'.format(numSeqs, assemblySize))
    else:
        status('Megahit assembly output missing -- check megahit logfile.')

    if not args.pipe:
        status('Your next command might be:\n\tAAFTF vecscreen -i {:} -c {:}\n'.format(finalOut, args.cpus))


def run(parser,args):
    if args.method  == "spades":
        run_spades(parser,args)
    elif args.method == "dipspades":
        run_dipspades(parser,args)
    elif args.method  == "megahit":
        run_megahit(parser,args)
    else:
        status("Unknow assembler method {}".format(args.method))
