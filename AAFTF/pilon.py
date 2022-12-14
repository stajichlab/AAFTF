"""Runs the pilon illumina-based polishing of assembly.

This takes care of running multiple rounds and updating the assembly
through these versions and removing temporary BAM alignment files.
"""

import os
import shutil
import subprocess
import sys
import uuid

from AAFTF.utility import SafeRemove, line_count, printCMD, status


def run(parser, args):
    """Runs pilon in multiple rounds provided with illumina reads and contig assembly FastA file."""
    # find reads for pilon
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)

    if not forReads:
        status('Unable to located FASTQ raw reads, ' +
               'pass via -l,--left and/or -r,--right')
        sys.exit(1)

    custom_workdir = 1
    if not args.workdir:
        custom_workdir = 0
        args.workdir = 'aaftf-pilon_'+str(uuid.uuid4())[:8]
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    bamthreads = 4
    if args.cpus < 4:
        bamthreads = args.cpus

    DEVNULL = open(os.devnull, 'w')
    for i in range(1, args.iterations+1):
        status(f'Starting Pilon polishing iteration {i}')
        correctedFasta = 'pilon'+str(i)+'.fasta'
        if i == 1:  # first loop
            initialFasta = args.infile
            shutil.copyfile(args.infile,
                            os.path.join(args.workdir,
                                         os.path.basename(args.infile)))
        else:
            initialFasta = os.path.join(args.workdir,
                                        'pilon'+str(i-1)+'.fasta')

        pilonBAM = os.path.basename(initialFasta)+'.bwa.bam'
        if not os.path.isfile(os.path.join(args.workdir, pilonBAM)):
            bwa_index = ['bwa', 'index', os.path.basename(initialFasta)]
            printCMD(bwa_index)
            subprocess.run(bwa_index, cwd=args.workdir, stderr=DEVNULL)
            bwa_cmd = ['bwa', 'mem', '-t', str(args.cpus),
                       os.path.basename(initialFasta), forReads]
            if revReads:
                bwa_cmd.append(revReads)

            # run BWA and pipe to samtools sort
            printCMD(bwa_cmd)
            p1 = subprocess.Popen(bwa_cmd, cwd=args.workdir,
                                  stdout=subprocess.PIPE, stderr=DEVNULL)
            p2 = subprocess.Popen(['samtools', 'sort',
                                   '-@', str(bamthreads), '-o', pilonBAM, '-'],
                                  cwd=args.workdir, stdout=subprocess.PIPE,
                                  stderr=DEVNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()

            # BAM file needs to be indexed for Pilon
            subprocess.run(['samtools', 'index', pilonBAM], cwd=args.workdir)

        # run Pilon
        pilon_cmd = ['pilon', '--genome', os.path.basename(initialFasta),
                     '--frags', pilonBAM,
                     f'-Xmx{args.memory}g',
                     '--output', correctedFasta.split('.fasta')[0],
                     '--threads', str(args.cpus),
                     '--changes']
        pilon_log = 'pilon'+str(i)+'.log'
        printCMD(pilon_cmd)
        with open(os.path.join(args.workdir, pilon_log), 'w') as logfile:
            subprocess.run(pilon_cmd, cwd=args.workdir, stderr=logfile,
                           stdout=logfile)
        n_chg = line_count(os.path.join(args.workdir,
                                        'pilon'+str(i)+'.changes'))

        status(f'Found {n_chg:,} changes in Pilon iteration {i}')

        # clean-up as we iterate to prevent tmp directory from blowing up
        dirty = [initialFasta+'.sa', initialFasta+'.amb', initialFasta+'.ann',
                 initialFasta+'.pac', initialFasta+'.bwt',
                 os.path.join(args.workdir, pilonBAM),
                 os.path.join(args.workdir, pilonBAM+'.bai')]
        for f in dirty:
            if i == 1:
                if os.path.isfile(os.path.join(args.workdir, f)):
                    os.remove(os.path.join(args.workdir, f))
            else:
                if os.path.isfile(f):
                    os.remove(f)

    # copy last iteration to output
    if args.outfile:
        polishedFasta = args.outfile
    else:
        fbasename = os.path.basename(args.infile).split('.f')[0]
        polishedFasta = f"{fbasename}.pilon.fasta"

    shutil.copyfile(os.path.join(args.workdir,
                                 'pilon'+str(args.iterations)+'.fasta'),
                    polishedFasta)

    status(f'AAFTF pilon completed {args.iterations} iterations.')
    status(f'Pilon polished assembly: {polishedFasta}')
    if '_' in polishedFasta:
        nextOut = polishedFasta.split('_')[0]+'.final.fasta'
    elif '.' in polishedFasta:
        nextOut = polishedFasta.split('.')[0]+'.final.fasta'
    else:
        nextOut = polishedFasta+'.final.fasta'

    if not args.debug and not custom_workdir:
        SafeRemove(args.workdir)

    if not args.pipe:
        status('Your next command might be:\n' +
               f'\tAAFTF sort -i {polishedFasta} -o {nextOut}\n')
