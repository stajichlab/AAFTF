"""Module to run a genome assembly using defaults for Fungi.

This uses SPAdes by default but additional tools like megahit are
supported and can be added. There is some access to updating
parameters but this entire package is intended to be a general
solution for draft Illumina genome processing en masse.
"""
import os
import re
import shutil
import subprocess
import sys
import uuid

from AAFTF.utility import fastastats, printCMD, status


def run_spades(parser, args):
    """Run SPAdes assembhler."""
    if not args.workdir:
        args.workdir = 'spades_'+str(uuid.uuid4())[:8]

    runcmd = ['spades.py', '--threads', str(args.cpus),
              '--mem', args.memory, '-o', args.workdir]

    if args.isolate:
        runcmd.extend(['--isolate'])
    elif args.careful:
        runcmd.extend(['--careful'])

    if args.assembler_args:
        runcmd.extend(args.assembler_args)

    if '--meta' not in runcmd:
        runcmd.extend(['--cov-cutoff', 'auto'])

    if args.tmpdir:
        runcmd.extend(['--tmp-dir', args.tmpdir])

    # find reads -- use --left/right or look for cleaned in tmpdir
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        status('Unable to located FASTQ raw reads, provide --left')
        sys.exit(1)

    if not revReads:
        runcmd .extend(['--s1', forReads])
        if args.merged:
            runcmd.extend(['--s2', args.merged])
    else:
        runcmd.extend(['--pe1-1', forReads, '--pe1-2', revReads])
        if args.merged:
            runcmd.extend(['--s1', args.merged])

    # this basically overrides everything above and only runs --restart-from option
    if os.path.isdir(args.workdir):
        runcmd = ['spades.py', '-o', args.workdir,
                  '--threads', str(args.cpus),
                  '--mem', args.memory,
                  '--restart-from last']

    # now run the spades job
    status('Assembling FASTQ data using Spades')
    printCMD(runcmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(runcmd)
    else:
        subprocess.run(runcmd, stdout=DEVNULL, stderr=DEVNULL)

    # pull out assembly
    if args.out:
        finalOut = args.out
    else:
        prefix = os.basename(forReads)
        m = re.search(r'(\S+)\.(fastq|fq)(\.\S+)?', prefix)
        if m:
            prefix = m.group(1)
        finalOut = prefix+'.spades.fasta'

    if os.path.isfile(os.path.join(args.workdir, 'scaffolds.fasta')):
        shutil.copyfile(os.path.join(args.workdir, 'scaffolds.fasta'), finalOut)
        status(f'Spades assembly finished: {finalOut}')
        numSeqs, assemblySize = fastastats(finalOut)
        status(f'Assembly is {numSeqs:,} scaffolds and {assemblySize:,} bp')
    else:
        status('Spades assembly output missing -- check Spades logfile.')

    if not args.pipe:
        status(f'Your next command might be:\n\tAAFTF vecscreen -i {finalOut} -c {args.cpus}\n')


def run_dipspades(parser, args):
    """Run dipSPAdes for diploid assembly support, only on older version of SPAdes."""
    if not args.workdir:
        args.workdir = 'dipspades_'+str(os.getpid())

    runcmd = ['dipspades.py', '--threads', str(args.cpus), '--cov-cutoff',
              'auto', '--mem', args.memory, '-o',
              args.workdir]

    if args.assembler_args:
        runcmd.extend(args.assembler_args)

    if args.haplocontigs:
        runcmd.extend(['--hap', args.haplocontigs])

    if args.tmpdir:
        runcmd.extend(['--tmp-dir', args.tmpdir])

    # find reads -- use --left/right or look for cleaned in tmpdir
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        status('Unable to located FASTQ raw reads, provide --left')
        sys.exit(1)

    if not revReads:
        runcmd.extend(['-s', forReads])
    else:
        runcmd.extend(['--pe1-1', forReads, '--pe1-2', revReads])
        if args.merged:
            runcmd.extend(['-s', args.merged])

    # this basically overrides everything above and only runs --restart-from option
    if os.path.isdir(args.workdir):
        runcmd = ['dipspades.py', '-o', args.workdir,
                  '--continue']

    # now run the spades job
    status('Assembling FASTQ data using Spades')

    printCMD(runcmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(runcmd)
    else:
        subprocess.run(runcmd, stdout=DEVNULL, stderr=DEVNULL)

    # pull out assembly file
    if args.out:
        finalOut = args.out
    else:
        prefix = os.basename(forReads)
        m = re.search(r'(\S+)\.(fastq|fq)(\.\S+)?', prefix)
        if m:
            prefix = m.group(1)
        finalOut = prefix+'.dipspades.fasta'

    if os.path.isfile(os.path.join(args.workdir, 'consensus_contigs.fasta')):
        shutil.copyfile(os.path.join(args.workdir, 'consensus_contigs.fasta'), finalOut)
        shutil.copyfile(os.path.join(args.workdir, 'dipspades', 'paired_consensus_contigs.fasta'),
                        prefix+".dipspades_consensus_paired.fasta")
        shutil.copyfile(os.path.join(args.workdir, 'dipspades', 'paired_consensus_contigs.fasta'),
                        prefix+".dipspades_consensus_unpaired.fasta")
        status(f'Dipspades assembly finished: {finalOut}')
        status('Dipspades assembly copied over: {:}'.format(prefix+".dipspades_consensus_unpaired.fasta"),
               prefix+".dipspades_consensus_paired.fasta")
        numSeqs, assemblySize = fastastats(finalOut)
        status(f'Assembly is {numSeqs:,} scaffolds and {assemblySize:,} bp')
    else:
        status('Spades assembly output missing -- '
               'check Dipspades logfile in {:}.'.format(os.path.join(args.workdir,
                                                                     'dipspades',
                                                                     'dipspades.log')))

    if not args.pipe:
        status(f'Your next command might be:\n\tAAFTF vecscreen -i {finalOut} -c {args.cpus}\n')


def run_megahit(parser, args):
    """Run megahit assembler. This is faster but maybe less accurate."""
    if not args.workdir:
        args.workdir = 'megahit_'+str(os.getpid())

    runcmd = ['megahit', '-t', str(args.cpus),
              '-o', args.workdir]

    if args.assembler_args:
        runcmd.extend(args.assembler_args)

    if args.memory:
        runcmd.extend(['--memory', args.memory])

    if args.tmpdir:
        runcmd.extend(['--tmp-dir', args.tmpdir])

    # find reads -- use --left/right or look for cleaned in tmpdir
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        status('Unable to located FASTQ raw reads, provide --left')
        sys.exit(1)

    if not revReads:
        runcmd.extend(['-r', forReads])
    else:
        runcmd.extend(['-1', forReads, '-2', revReads])

    if os.path.isdir(args.workdir):
        status(f"Cannot re-run with existing folder {args.workdir}")

    # now run the spades job
    status('Assembling FASTQ data using megahit')
    printCMD(runcmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(runcmd)
    else:
        subprocess.run(runcmd, stdout=DEVNULL, stderr=DEVNULL)
    # pull out assembly
    if args.out:
        finalOut = args.out
    else:
        prefix = os.basename(forReads)
        m = re.search(r'(\S+)\.(fastq|fq)(\.\S+)?', prefix)
        if m:
            prefix = m.group(1)
        finalOut = prefix+'.megahit.fasta'

    if os.path.isfile(os.path.join(args.workdir, 'final.contigs.fa')):
        shutil.copyfile(os.path.join(args.workdir, 'final.contigs.fa'), finalOut)
        status(f'Megahit assembly finished: {finalOut}')
        numSeqs, assemblySize = fastastats(finalOut)
        status(f'Assembly is {numSeqs:,} scaffolds and {assemblySize:,} bp')
    else:
        status('Megahit assembly output missing -- check megahit logfile.')

    if not args.pipe:
        status(f'Your next command might be:\n\tAAFTF vecscreen -i {finalOut} -c {args.cpus}\n')


def run_unicycler(parser, args):
    """Run Unicycler assembhler."""
    if not args.workdir:
        args.workdir = 'unicycler_'+str(uuid.uuid4())[:8]

    runcmd = ['unicycler', '--threads', str(args.cpus),
              '-o', args.workdir]

    # if args.memory:
    #    runcmd.extend(['--spades_options', f'-m {args.memory}'])

    # find reads -- use --left/right or look for cleaned in tmpdir
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        status('Unable to located FASTQ raw reads, provide --left')
        sys.exit(1)

    if args.longreads:
        runcmd.extend(['--long', args.longreads])

    if not revReads:
        runcmd.extend(['--unpaired', forReads])
    elif args.merged:
        runcmd.extend(['--unpaired', args.merged])
    else:
        runcmd.extend(['--short1', forReads, '--short2', revReads])
        if args.merged:
            runcmd.extend(['--unpaired', args.merged])

    # not supporting restarting a run
    # this basically overrides everything above and only runs --restart-from option
    #    if os.path.isdir(args.workdir):
    #    runcmd = ['unicycler', '-o', args.workdir,
    #            '--threads', str(args.cpus),
    #            '--mem', args.memory,
    #            '--restart-from last']

    # now run the spades job
    status('Assembling FASTQ data using Unicycler')
    printCMD(runcmd)
    DEVNULL = open(os.devnull, 'w')
    if args.debug:
        subprocess.run(runcmd)
    else:
        subprocess.run(runcmd, stdout=DEVNULL, stderr=DEVNULL)

    # pull out assembly
    if args.out:
        finalOut = args.out
    else:
        prefix = os.basename(forReads)
        m = re.search(r'(\S+)\.(fastq|fq)(\.\S+)?', prefix)
        if m:
            prefix = m.group(1)
        finalOut = prefix+'.unicycler.fasta'

    if os.path.isfile(os.path.join(args.workdir, 'assembly.fasta')):
        shutil.copyfile(os.path.join(args.workdir, 'assembly.fasta'), finalOut)
        status(f'Unicycler assembly finished: {finalOut}')
        numSeqs, assemblySize = fastastats(finalOut)
        status(f'Assembly is {numSeqs:,} scaffolds and {assemblySize:,} bp')
    else:
        status('Unicycler assembly output missing -- check Unicycler logfile.')

    if not args.pipe:
        status(f'Your next command might be:\n\tAAFTF vecscreen -i {finalOut} -c {args.cpus}\n')


def run(parser, args):
    """General run command for this subcommand module where parameters are consumed."""
    if args.method == "spades":
        run_spades(parser, args)
    elif args.method == "dipspades":
        run_dipspades(parser, args)
    elif args.method == "megahit":
        run_megahit(parser, args)
    elif args.method == "masurca":
        status("Masurca assembly is not yet implemented in AAFTF")
    elif args.method == "nextdenovo":
        status("NextDenovo assembly is not yet implemented in AAFTF")
    elif args.method == "unicycler":
        run_unicycler(parser, args)
    else:
        status(f"Unknown assembler method {args.method}")
