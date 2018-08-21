import sys, os, shutil, subprocess

#logging
import logging
logger = logging.getLogger('AAFTF')

from AAFTF.utility import line_count

def run(parser,args):

    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)
    
    #find reads for pilon
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
        logger.error('Unable to located FASTQ raw reads')
        sys.exit(1)
        
    if not args.prefix:
        if '_' in os.path.basename(forReads):
            args.prefix = os.path.basename(forReads).split('_')[0]
        elif '.' in os.path.basename(forReads):
            args.prefix = os.path.basename(forReads).split('.')[0]
        else:
            args.prefix = os.path.basename(forReads)
    
    DEVNULL = open(os.devnull, 'w')
    for i in range(1, args.iterations+1):
        logger.info(' Starting Pilon polishing iteration {:}'.format(i))
        correctedFasta = 'pilon'+str(i)+'.fasta'
        if i == 1: #first loop
            initialFasta = args.infile
            shutil.copyfile(args.infile, os.path.join(args.workdir, os.path.basename(args.infile)))
        else:
            initialFasta = os.path.join(args.workdir, 'pilon'+str(i-1)+'.fasta')
        pilonBAM = os.path.basename(initialFasta)+'.bwa.bam'
        if not os.path.isfile(os.path.join(args.workdir, pilonBAM)):
            bwa_index = ['bwa', 'index', os.path.basename(initialFasta)]
            logger.info('CMD: {:}'.format(' '.join(bwa_index)))
            subprocess.run(bwa_index, cwd=args.workdir, stderr=DEVNULL)
            bwa_cmd = ['bwa', 'mem', '-t', str(args.cpus), os.path.basename(initialFasta), forReads]
            if revReads:
                bwa_cmd.append(revReads)
    
            #run BWA and pipe to samtools sort
            logger.info(' CMD: {:}'.format(' '.join(bwa_cmd)))
            p1 = subprocess.Popen(bwa_cmd, cwd=args.workdir, stdout=subprocess.PIPE, stderr=DEVNULL)
            p2 = subprocess.Popen(['samtools', 'sort', '-@', str(args.cpus),'-o', pilonBAM, '-'], 
                        cwd=args.workdir, stdout=subprocess.PIPE, stderr=DEVNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()
            
            #BAM file needs to be indexed for Pilon
            subprocess.run(['samtools', 'index', pilonBAM], cwd=args.workdir)
        
        #run Pilon
        pilon_cmd = ['pilon', '--genome', os.path.basename(initialFasta), '--frags', pilonBAM, 
                     '--output', correctedFasta.split('.fasta')[0], '--threads', str(args.cpus), 
                     '--changes']
        pilon_log = 'pilon'+str(i)+'.log'
        logger.info(' CMD: {:}'.format(' '.join(pilon_cmd)))
        with open(os.path.join(args.workdir, pilon_log), 'w') as logfile:
            subprocess.run(pilon_cmd, cwd=args.workdir, stderr=logfile, stdout=logfile)
        num_changes = line_count(os.path.join(args.workdir, 'pilon'+str(i)+'.changes'))
        logger.info(' found {:,} changes in Pilon iteration {:}'.format(num_changes, i))
        
        #clean-up as we iterate to prevent tmp directory from blowing up
        dirty = [initialFasta+'.sa', initialFasta+'.amb', initialFasta+'.ann',
                 initialFasta+'.pac', initialFasta+'.bwt', os.path.join(args.workdir, pilonBAM), 
                 os.path.join(args.workdir, pilonBAM+'.bai')]
        for f in dirty:
            if i == 1:
                if os.path.isfile(os.path.join(args.workdir, f)):
                    os.remove(os.path.join(args.workdir, f))
            else:
                if os.path.isfile(f):
                    os.remove(f)
    
    #copy last iteration to output
    if args.outfile:
        polishedFasta = args.outfile
    else:
        polishedFasta = os.path.basename(args.infile).split('.f')[0]+'.pilon.fasta'
    shutil.copyfile(os.path.join(args.workdir, 'pilon'+str(args.iterations)+'.fasta'), polishedFasta)
    logger.info(' AAFTF pilon completed {:} iterations.'.format(args.iterations))
    logger.info(' Pilon polished assembly: {:}'.format(polishedFasta))
    if '_' in polishedFasta:
        nextOut = polishedFasta.split('_')[0]+'.final.fasta'
    elif '.' in polishedFasta:
        nextOut = polishedFasta.split('.')[0]+'.final.fasta'
    else:
        nextOut = polishedFasta+'.final.fasta'
    logger.info('Your next command might be:\n\tAAFTF sort -i {:} -o {:}\n'.format(polishedFasta, nextOut))      
        
