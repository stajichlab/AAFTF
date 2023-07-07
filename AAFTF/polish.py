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


def run(parser, args):  # noqa: C901
    """Execute polishing step in multiple rounds provided with illumina reads and contig assembly FastA file."""
    # find reads for polishing

    polishMethod = args.method
    memperthread = int(args.memory / args.cpus)
    if memperthread == 0:
        memperthread = "500K"
    else:
        memperthread = f"{memperthread}G"

    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)

    longreads = None
    if args.longreads:
        # racon or nextpolish can use these
        longreads = os.path.abspath(args.longreads)
    elif polishMethod == "racon":
        status(f"calling {polishMethod} without long reads (pacbio/ONT)")
        sys.exit(1)

    if args.method == "racon" and not longreads:
        status('Unable to located long read FASTQ raw reads, ' +
               'pass via -lr or --longreads')
        sys.exit(1)
    if not forReads and args.method != "racon":
        status('Unable to located FASTQ raw reads, ' +
               'pass via -l,--left and/or -r,--right')
        sys.exit(1)

    custom_workdir = 1
    if not args.workdir:
        custom_workdir = 0
        args.workdir = 'aaftf-polish_'+str(uuid.uuid4())[:8]
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    # Output file
    polishedFasta = None
    if args.outfile:
        polishedFasta = args.outfile
    else:
        fbasename = os.path.basename(args.infile).split('.f')[0]
        polishedFasta = f"{fbasename}.polished.fasta"

    method = args.method.lower()
    nextPolishExe = None
    polish_log = f'{method}.log'
    if method == "pilon" or method == "nextpolish":
        for i in range(1, args.iterations+1):
            status(f'Starting {method} polishing iteration {i}')
            correctedBase = f'polished{i}'
            if i == 1:  # first loop
                initialFasta = args.infile
                initialFasta = os.path.join(args.workdir,
                                            os.path.basename(args.infile))
                shutil.copyfile(args.infile, initialFasta)
            else:
                initialFasta = os.path.join(args.workdir,
                                            'polished'+str(i-1)+'.fasta')
            BAMfile = make_bwa_bam(initialFasta, forReads, revReads,
                                   args.workdir, args.cpus, memperthread)
            run_cmd = []
            polish_log = None
            dirty = []
            if method == "pilon":
                # run Pilon
                run_cmd = ['pilon', '--genome', os.path.basename(initialFasta),
                           '--frags', BAMfile,
                           f'-Xmx{args.memory}g',
                           '--output', correctedBase,
                           '--threads', str(args.cpus),
                           '--changes']
                if args.diploid or args.ploidy == 2:
                    run_cmd.append('--diploid')

                polish_log = 'pilon_'+str(i)+'.log'

                printCMD(run_cmd)
                with open(os.path.join(args.workdir, polish_log), 'w') as logfile:
                    subprocess.run(run_cmd, cwd=args.workdir,
                                   stderr=logfile, stdout=logfile)
                n_chg = line_count(os.path.join(args.workdir, correctedBase + '.changes'))

                status(f'Found {n_chg:,} changes in Pilon iteration {i}')
            elif method == "nextpolish":
                if not nextPolishExe:
                    nextPolishmain = shutil.which('nextPolish')
                    print(nextPolishmain)
                    nextPolishExe = os.path.join(os.path.dirname(os.path.dirname(nextPolishmain)),
                                                 'share', 'nextpolish-1.4.1', 'lib', 'nextpolish1.py')
                    if not os.path.exists(nextPolishExe):
                        nextPolishExe = os.path.join(os.path.dirname(nextPolishmain),
                                                     'lib', 'nextpolish1.py')
                if not nextPolishExe or not os.path.exists(nextPolishExe):
                    status("Cannot find nextPolish python script")
                    return -1
                print(f'np is {nextPolishExe}')
                # initialFasta should have already been copied to the working dir
                # or is carryforward from last iteration
                run_cmd = ['samtools', 'faidx', os.path.basename(initialFasta)]
                subprocess.run(run_cmd, cwd=args.workdir)
                tempoutfasta = f'temp_{correctedBase}.fasta'
                run_cmd = ['python', nextPolishExe,
                           '-g', os.path.basename(initialFasta),
                           '-t', '1',
                           '-s', BAMfile,
                           '-p', str(args.cpus),
                           '-o', tempoutfasta]
                if args.diploid or args.ploidy == 2:
                    run_cmd.extend(['-ploidy', '2'])
                elif args.ploidy:
                    run_cmd.extend(['-ploidy', str(args.ploidy)])
                polish_log = 'nextpolish_t1_'+str(i)+'.log'
                with open(os.path.join(args.workdir, polish_log), 'w') as logfile:
                    printCMD(run_cmd)
                    subprocess.run(run_cmd, cwd=args.workdir, stderr=logfile, stdout=logfile)
                logfile.close()
                run_cmd = ['samtools', 'faidx', tempoutfasta]
                # make second BAM file for second task of nextPolish
                BAMfile = make_bwa_bam(tempoutfasta, forReads, revReads,
                                       args.workdir, args.cpus, memperthread)

                run_cmd = ['python', nextPolishExe,
                           '-g', tempoutfasta,
                           '-t', '2',
                           '-debug',
                           '-p', str(args.cpus),
                           '-s', BAMfile,
                           '-o', correctedBase + ".fasta"]
                if args.diploid or args.ploidy == 2:
                    run_cmd.extend(['-ploidy', '2'])
                elif args.ploidy:
                    run_cmd.extend(['-ploidy', str(args.ploidy)])
                polish_log = 'nextpolish_t2_'+str(i)+'.log'
                with open(os.path.join(args.workdir, polish_log), 'w') as logfile:
                    printCMD(run_cmd)
                    subprocess.run(run_cmd, cwd=args.workdir, stderr=logfile, stdout=logfile)
                dirty.append(tempoutfasta)

            # clean-up as we iterate to prevent tmp directory from blowing up
            dirty.extend([initialFasta+'.sa', initialFasta+'.amb',
                          initialFasta+'.ann', initialFasta+'.pac',
                          initialFasta+'.bwt',
                          os.path.join(args.workdir, BAMfile),
                          os.path.join(args.workdir, BAMfile+'.bai')])
            for f in dirty:
                if i == 1:
                    if os.path.isfile(os.path.join(args.workdir, f)):
                        os.remove(os.path.join(args.workdir, f))
                else:
                    if os.path.isfile(f):
                        os.remove(f)

        shutil.copyfile(os.path.join(args.workdir, 'polished' + str(args.iterations)+'.fasta'),
                        polishedFasta)

        status(f'AAFTF polish completed {args.iterations} iterations.')
        status(f'{method} polished assembly: {polishedFasta}')
        if '_' in polishedFasta:
            nextOut = polishedFasta.split('_')[0]+'.final.fasta'
        elif '.' in polishedFasta:
            nextOut = polishedFasta.split('.')[0]+'.final.fasta'
        else:
            nextOut = polishedFasta+'.final.fasta'
    elif args.method.lower() == "polca" or args.method.lower() == "masurca":
        initialFasta = os.path.basename(args.infile)
        shutil.copyfile(args.infile, os.path.join(args.workdir, initialFasta))
        polca_cmd = [args.polca, '-a', initialFasta,
                     '-r', f'{forReads} {revReads}',
                     '-t', str(args.cpus), '-m', memperthread]
        printCMD(polca_cmd)
        # run the polca polishing
        with open(os.path.join(args.workdir, polish_log), 'w') as logfile:
            subprocess.run(polca_cmd, cwd=args.workdir,
                           stderr=logfile,
                           stdout=logfile)
        shutil.copyfile(os.path.join(args.workdir, f'{initialFasta}.PolcaCorrected.fa'),
                        polishedFasta)
        shutil.copyfile(os.path.join(args.workdir, f'{initialFasta}.vcf'),
                        f'{polishedFasta}.vcf')
        shutil.copyfile(os.path.join(args.workdir, f'{initialFasta}.report'),
                        f'{polishedFasta}.polca_report.txt')
        status('AAFTF polish completed.')
        status(f'{method} polished assembly: {polishedFasta}')

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


def make_bwa_bam(inFasta, forReads, revReads,
                 workdir, cpus, memperthread):
    """Run BAM file generation from short reads on current assembly file to enable polishing."""
    ASMname = os.path.basename(inFasta)
    ASMpref = os.path.splitext(ASMname)[0]
    BAM = ASMpref + '.bwa.bam'
    if os.path.exists(os.path.join(workdir, BAM)):
        return BAM
    tempfiles = [f'{ASMpref}.fixmate.bam', f'{ASMpref}.markdup.bam']
    DEVNULL = open(os.devnull, 'w')
    bamthreads = 4
    if cpus < 4:
        bamthreads = cpus
    if not os.path.isfile(os.path.join(workdir, BAM)):
        bwa_index = ['bwa', 'index', ASMname]
        printCMD(bwa_index)
        subprocess.run(bwa_index, cwd=workdir, stderr=DEVNULL)
        bwa_cmd = ['bwa', 'mem', '-t', str(cpus), ASMname, forReads]
        if revReads:
            bwa_cmd.append(revReads)

        # run BWA and pipe to samtools sort
        printCMD(bwa_cmd)
        p1 = subprocess.Popen(bwa_cmd, cwd=workdir,
                              stdout=subprocess.PIPE, stderr=DEVNULL)

        samtools_cmd = ['samtools', 'fixmate', '-O', 'bam,level=1', '-m', '-', tempfiles[0]]
        printCMD(samtools_cmd)
        p2 = subprocess.Popen(samtools_cmd, stdin=p1.stdout, cwd=workdir)
        p1.stdout.close()
        p2.communicate()

        # level 1
        samtools_cmd = ['samtools', 'sort', '-l', '1',
                        '-@', str(bamthreads), '-T', ASMpref,
                        '-m', memperthread, tempfiles[0]]
        printCMD(samtools_cmd)
        p3 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, cwd=workdir)

        samtools_cmd = ['samtools', 'markdup', '-@',  str(bamthreads), '-', tempfiles[1]]
        printCMD(samtools_cmd)
        p4 = subprocess.Popen(samtools_cmd, stdin=p3.stdout, cwd=workdir)
        p3.stdout.close()
        p4.communicate()

        # keep only paired reads
        samtools_cmd = ['samtools', 'view', '-O', 'bam,level=1', '-f', '0x2',
                        '-@',  str(bamthreads), '-o', BAM, tempfiles[1]]
        printCMD(samtools_cmd)
        subprocess.run(samtools_cmd, cwd=workdir)

        # BAM file needs to be indexed
        samtools_cmd = ['samtools', 'index', '-@', str(cpus), BAM]
        printCMD(samtools_cmd)
        subprocess.run(samtools_cmd, cwd=workdir)

        for tfile in tempfiles:
            os.remove(os.path.join(workdir, tfile))
    return BAM
