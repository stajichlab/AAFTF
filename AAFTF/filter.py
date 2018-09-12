import sys, os, shutil, gzip, subprocess
import urllib.request

# this runs rountines to remove sequence reads
# which match contaminant databases and sources
# including by default PhiX and others specified by
# the user

# it will download the libraries from GenBank using
# accession numbers as well as some hard coded links
# to PhiX - see resouces.py for these defaults

#logging
import logging
logger = logging.getLogger('AAFTF')

from AAFTF.resources import Contaminant_Accessions
from AAFTF.resources import SeqDBs
from AAFTF.resources import DB_Links
from AAFTF.utility import bam_read_count
from AAFTF.utility import countfastq

def run(parser,args):
    
    if args.workdir == 'working_AAFTF' and not args.prefix and not args.left and args.cpus == 1: 
        logger.error(' Please provide either -w,--workdir, -p,--prefix, or --left reads.')
        sys.exit(1)

    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    #parse database locations
    DB = None
    if not args.AAFTF_DB:
        try:
            DB = os.environ["AAFTF_DB"]
        except KeyError:
            if args.AAFTF_DB:
                DB = args.AAFTF_DB
            else:
                pass
    else:
        DB = args.AAFTF_DB
            
    earliest_file_age = -1
    contam_filenames = []
    # db of contaminant (PhiX)
    for url in Contaminant_Accessions.values():
        acc = os.path.basename(url)
        if DB:
            acc_file = os.path.join(DB, acc)
        else:
            acc_file = os.path.join(args.workdir,acc)
        contam_filenames.append(acc_file)
        if not os.path.exists(acc_file):
            urllib.request.urlretrieve(url,acc_file)
        if ( earliest_file_age < 0 or
             earliest_file_age < os.path.getctime(acc_file) ):
            earliest_file_age = os.path.getctime(acc_file)

    # download univec too
    url = DB_Links['UniVec']
    acc = os.path.basename(DB_Links['UniVec'])
    if DB:
        acc_file = os.path.join(DB, acc)
    else:
        acc_file = os.path.join(args.workdir,acc)
    contam_filenames.append(acc_file)
    if not os.path.exists(acc_file):
        urllib.request.urlretrieve(url,acc_file)
        if ( earliest_file_age < 0 or
             earliest_file_age < os.path.getctime(acc_file) ):
            earliest_file_age = os.path.getctime(acc_file)
    
    if args.screen_accessions:
        for acc in args.screen_accessions:
            if DB:
                acc_file = os.path.join(DB, acc+".fna")
                if not os.path.exists(acc_file):
                    acc_file = os.path.join(args.workdir,acc+".fna")
            else:
                acc_file = os.path.join(args.workdir,acc+".fna")
            contam_filenames.append(acc_file)
            if not os.path.exists(acc_file):
                url = SeqDBs['nucleotide'] % (acc)
                urllib.request.urlretrieve(url,acc_file)
            if ( earliest_file_age < 0 or
                 earliest_file_age < os.path.getctime(acc_file) ):
                earliest_file_age = os.path.getctime(acc_file)

    if args.screen_urls:
        for url in args.screen_urls:
            url_file = os.path.join(args.workdir,os.path.basename(url))
            contam_filenames.append(url_file)
            if not os.path.exists(url_file):
                urllib.request.urlretrieve(url,url_file)
            if ( earliest_file_age < 0 or
                 earliest_file_age < os.path.getctime(url_file) ):
                earliest_file_age = os.path.getctime(url_file)

    # concat vector db
    logger.info(' Generating combined contamination database:\n{:}'.format('\n'.join(contam_filenames)))
    contamdb = os.path.join(args.workdir,'contamdb')
    if ( not os.path.exists(contamdb) or
         ( os.path.getctime(contamdb) < earliest_file_age)):
         with open(contamdb, 'wb') as wfd:
             for fname in contam_filenames:
                 with open(fname,'rb') as fd: # reasonably fast copy for append
                     shutil.copyfileobj(fd, wfd)
    
    #find reads
    forReads, revReads = (None,)*2
    if args.left:
        forReads = os.path.abspath(args.left)
    if args.right:
        revReads = os.path.abspath(args.right)
    if not forReads:
        if not args.prefix:
            logger.info("Must provide prefix if no left/right file pairs")
        else:
            for file in os.listdir(args.workdir):
                if file.startswith(args.prefix) and '_1P.fastq' in file:
                    forReads = os.path.abspath(os.path.join(args.workdir, file))
                if file.startswith(args.prefix) and '_2P.fastq' in file:
                    revReads = os.path.abspath(os.path.join(args.workdir, file))
    if not forReads:
        logger.error(' Unable to located FASTQ reads')
        sys.exit(1)
    
    # seems like this needs to be stripping trailing extension?
    if not args.prefix:
        if '_' in os.path.basename(forReads):
            args.prefix = os.path.basename(forReads).split('_')[0]
        elif '.' in os.path.basename(forReads):
            args.prefix = os.path.basename(forReads).split('.')[0]
        else:
            args.prefix = os.path.basename(forReads)
        
    #logger.info('Loading {:,} FASTQ reads'.format(countfastq(forReads)))
    DEVNULL = open(os.devnull, 'w')
    alignBAM = os.path.join(args.workdir, args.prefix+'_contam_db.bam')
    clean_reads = os.path.join(args.workdir, args.prefix + "_cleaned")
    
    if args.aligner == 'bowtie2':  
        if not os.path.isfile(alignBAM):
            logger.info(' Aligning reads to contamination database using bowtie2')
            if ( not os.path.exists(contamdb + ".1.bt2") or
                 os.path.getctime(contamdb + ".1.bt2") < 
                 os.path.getctime(contamdb)):
                # (re)build index if no index or index is older than 
                # the db
                bowtie_index = ['bowtie2-build', contamdb, contamdb]
                logger.info('CMD: {:}'.format(' '.join(bowtie_index)))
                subprocess.run(bowtie_index, stderr=DEVNULL, stdout=DEVNULL)

            bowtie_cmd = ['bowtie2','-x', os.path.basename(contamdb),
                          '-p', str(args.cpus), '--very-sensitive']
            if forReads and revReads:
                bowtie_cmd = bowtie_cmd + ['-1', forReads, '-2', revReads]
            elif forReads:
                bowtie_cmd = bowtie_cmd + ['-U', forReads]
            
            #now run and write to BAM sorted
            logger.info('CMD: {:}'.format(' '.join(bowtie_cmd)))
            p1 = subprocess.Popen(bowtie_cmd, cwd=args.workdir, stdout=subprocess.PIPE, stderr=DEVNULL)
            p2 = subprocess.Popen(['samtools', 'sort', '-@', str(args.cpus),
                                   '-o', os.path.basename(alignBAM), '-'],
                                   cwd=args.workdir, stdout=subprocess.PIPE, 
                                   stderr=DEVNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()        
                
    elif args.aligner == 'bwa':
        if not os.path.isfile(alignBAM):
            logger.info(' Aligning reads to contamination database using BWA')
            if ( not os.path.exists(contamdb + ".amb") or
                 os.path.getctime(contamdb + ".amb") < 
                 os.path.getctime(contamdb)):
                bwa_index = ['bwa','index', contamdb]
                logger.info('CMD: {:}'.format(' '.join(bwa_index)))
                subprocess.run(bwa_index, stderr=DEVNULL, stdout=DEVNULL)
            
            bwa_cmd = ['bwa', 'mem', '-t', str(args.cpus), os.path.basename(contamdb), forReads]
            if revReads:
                bwa_cmd.append(revReads)
            
            #now run and write to BAM sorted
            logger.info('CMD: {:}'.format(' '.join(bwa_cmd)))
            p1 = subprocess.Popen(bwa_cmd, cwd=args.workdir, stdout=subprocess.PIPE, stderr=DEVNULL)
            p2 = subprocess.Popen(['samtools', 'sort', '-@', str(args.cpus),
                                   '-o', os.path.basename(alignBAM), '-'],
                                   cwd=args.workdir, stdout=subprocess.PIPE, 
                                   stderr=DEVNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()                
      
    elif args.aligner == 'minimap2':
        if not os.path.isfile(alignBAM):
            logger.info(' Aligning reads to contamination database using minimap2')
            
            minimap2_cmd = ['minimap2', '-ax', 'sr', '-t', str(args.cpus), os.path.basename(contamdb), forReads]
            if revReads:
                minimap2_cmd.append(revReads)
            
            #now run and write to BAM sorted
            logger.info('CMD: {:}'.format(' '.join(minimap2_cmd)))
            p1 = subprocess.Popen(minimap2_cmd, cwd=args.workdir, stdout=subprocess.PIPE, stderr=DEVNULL)
            p2 = subprocess.Popen(['samtools', 'sort', '-@', str(args.cpus),
                                   '-o', os.path.basename(alignBAM), '-'],
                                   cwd=args.workdir, stdout=subprocess.PIPE, 
                                   stderr=DEVNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()    
    else:
        logger.error(" Must specify bowtie2, bwa, or minimap2 for filtering")
    
    if os.path.isfile(alignBAM):
        #display mapping stats in terminal
        subprocess.run(['samtools', 'index', alignBAM])
        mapped, unmapped = bam_read_count(alignBAM)
        logger.info(' {:,} reads mapped to contamintation database'.format(mapped))
        logger.info(' {:,} reads unmapped and writing to file'.format(unmapped))
        #now output unmapped reads from bamfile
        #this needs to be -f 5 so unmapped-pairs
        if forReads and revReads:
            samtools_cmd = ['samtools', 'fastq', '-f', '12',
                            '-1', clean_reads+'_1.fastq.gz',
                            '-2', clean_reads+'_2.fastq.gz',
                            alignBAM]
        elif forReads:
            samtools_cmd = ['samtools', 'fastq', '-f', '4',
                            '-1', clean_reads+'_1.fastq.gz',
                            alignBAM]
        subprocess.run(samtools_cmd, stderr=DEVNULL)
        if not args.debug:
            os.remove(alignBAM)
            os.remove(alignBAM + ".bai")
        if revReads:
            logger.info(' Filtering complete:\n\tFor: {:}\n\tRev: {:}'.format(
                        clean_reads+'_1.fastq.gz',clean_reads+'_2.fastq.gz'))
        else:
            logger.info(' Filtering complete:\n\tSingle: {:}'.format(clean_reads+'_1.fastq.gz'))
        logger.info('Your next command might be:\n\tAAFTF assemble -w {:} -c {:} -o {:}\n'.format(
                    args.workdir, args.cpus, args.prefix+'.spades.fasta'))
            
