import sys, os, shutil, subprocess

#logging
import logging
logger = logging.getLogger('AAFTF')

from AAFTF.utility import line_count

def run(parser,args):
    
    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    
    #find reads for pilon
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        for file in os.listdir(args.tmpdir):
            if '_cleaned' in file and file.endswith('.fq.gz'):
                if '_1' in file:
                    forReads = os.path.abspath(os.path.join(args.tmpdir, file))
                if '_2' in file:
                    revReads = os.path.abspath(os.path.join(args.tmpdir, file))
    if not forReads:
        print('Unable to located FASTQ raw reads')
        sys.exit(1)
    
    DEVNULL = open(os.devnull, 'w')
    for i in range(1, args.iterations+1):
        logger.info(' Starting Pilon polishing iteration {:}'.format(i))
        correctedFasta = 'pilon'+str(i)+'.fasta'
        if i == 1: #first loop
            initialFasta = args.infile
            shutil.copyfile(args.infile, os.path.join(args.tmpdir, os.path.basename(args.infile)))
        else:
            initialFasta = os.path.join(args.tmpdir, 'pilon'+str(i-1)+'.fasta')
        pilonBAM = os.path.basename(initialFasta)+'.bwa.bam'
        if not os.path.isfile(os.path.join(args.tmpdir, pilonBAM)):
            bwa_index = ['bwa', 'index', os.path.basename(initialFasta)]
            print(' CMD: {:}'.format(' '.join(bwa_index)))
            subprocess.run(bwa_index, cwd=args.tmpdir, stderr=DEVNULL)
            bwa_cmd = ['bwa', 'mem', '-t', str(args.cpus), os.path.basename(initialFasta), forReads]
            if revReads:
                bwa_cmd.append(revReads)
    
            #run BWA and pipe to samtools sort
            print(' CMD: {:}'.format(' '.join(bwa_cmd)))
            p1 = subprocess.Popen(bwa_cmd, cwd=args.tmpdir, stdout=subprocess.PIPE, stderr=DEVNULL)
            p2 = subprocess.Popen(['samtools', 'sort', '-@', str(args.cpus),'-o', pilonBAM, '-'], 
            			cwd=args.tmpdir, stdout=subprocess.PIPE, stderr=DEVNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()
            
            #BAM file needs to be indexed for Pilon
            subprocess.run(['samtools', 'index', pilonBAM], cwd=args.tmpdir)
        
        #run Pilon
        pilon_cmd = ['pilon', '--genome', os.path.basename(initialFasta), '--frags', pilonBAM, 
                     '--output', correctedFasta.split('.fasta')[0], '--threads', str(args.cpus), 
                     '--changes']
        pilon_log = 'pilon'+str(i)+'.log'
        print(' CMD: {:}'.format(' '.join(pilon_cmd)))
        with open(os.path.join(args.tmpdir, pilon_log), 'w') as logfile:
            subprocess.run(pilon_cmd, cwd=args.tmpdir, stderr=logfile, stdout=logfile)
        num_changes = line_count(os.path.join(args.tmpdir, 'pilon'+str(i)+'.changes'))
        logger.info(' found {:,} changes in Pilon iteration {:}'.format(num_changes, i))
        
        #clean-up as we iterate to prevent tmp directory from blowing up
        dirty = [initialFasta+'.sa', initialFasta+'.amb', initialFasta+'.ann',
                 initialFasta+'.pac', initialFasta+'.bwt', os.path.join(args.tmpdir, pilonBAM), 
                 os.path.join(args.tmpdir, pilonBAM+'.bai')]
        for f in dirty:
            if i == 1:
                if os.path.isfile(os.path.join(args.tmpdir, f)):
                    os.remove(os.path.join(args.tmpdir, f))
            else:
                if os.path.isfile(f):
                    os.remove(f)
    
    #copy last iteration to output
    if args.outfile:
        polishedFasta = args.outfile
    else:
        polishedFasta = os.path.basename(args.infile).split('.f')[0]+'.pilon.fasta'
    shutil.copyfile(os.path.join(args.tmpdir, 'pilon'+str(args.iterations)+'.fasta'), polishedFasta)
    logger.info(' AAFTF pilon completed {:} iterations.'.format(args.iterations))
    logger.info(' Pilon polished assembly: {:}'.format(polishedFasta))
        
        