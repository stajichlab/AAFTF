import sys, os
import shutil

from subprocess import call, run, Popen, PIPE, STDOUT

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#logging

import logging
logger = logging.getLogger('AAFTF')


def run(parser,args):

    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)

    #find reads for blobplot
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
        
    # hard coded tmpfile
    assembly_working = 'assembly.fasta'
    megablast_working = 'megablast.out'
    blobBAM="remapped.bam"
    shutils.copy(args.input,os.path.join(args.tmpdir,assembly_working))
    # index
    bwa_index  = ['bwa','index',assembly_working]
    
    if args.debug:
        logger.debug(bwa_index)
        
    subprocess.run(bwa_index, cwd=args.tmpdir, stderr=DEVNULL)

    # get ready to map reads        
    # bwa is forcing these to be 2 pairs not interleaved
    # may need to detect/correct this if not
    left = os.path.join(args.indir,args.prefix + "_cleaned_1.fq.gz")
    right = os.path.join(args.indir,args.prefix + "_cleaned_2.fq.gz")
    #mapped reads to assembly using BWA
    bwa_args = ['bwa','mem',
               '-t', str(args.cpus),
                assembly_working, # assembly index base
                forReads]
    
    if revReads:
        bwa_args.append(revReads)
    
    #run BWA and pipe to samtools sort

    logger.info('CMD: {:}'.format(' '.join(bwa_cmd)))
    if os.path.exists(os.path.join(args.tmpdir,
                                   blobBLAM)):
        logger.info("BAM file %s already exists, not rerunning" %
                    os.path.join(args.tmpdir,blobBLAM))
    else:
        p1 = subprocess.Popen(bwa_args, cwd=args.tmpdir,
                              stdout=subprocess.PIPE, stderr=DEVNULL)
        p2 = subprocess.Popen(['samtools', 'sort', '--threads', str(args.cpus),
                               '-o', blobBAM, '-'], cwd=args.tmpdir,
                              stdout=subprocess.PIPE, stderr=DEVNULL,
                              stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()

        subprocess.run(['samtools', 'index', blobBAM], cwd=args.tmpdir)

    #resulting assembly was blasted against NCBI nt database
    blastn_cmd = ['blastn','-task', 'megablast',
                  '-query',assembly_working,
                  '-db', args.blastdb,
                  '-outfmt', "'6 qseqid staxids bitscore std sscinames sskingdoms stitle'",
                  '-culling_limit', '5',
                  '-num_threads', str(args.cpus),
                  '-evalue',args.evalue,
                  '-out', megablast_working]
    if args.debug:
        logger.debug(blastn_cmd)
    if os.path.exists(os.path.join(args.tmpdir,megablas_working)):
        logger.info("Megablast out %s already exists, not rerunning" % (os.path.join(args.tmpdir,megablas_working)))
    else:
        subprocess.run(blastn_cmd,cwd=args.tmpdir)

    blob_cmd = ['blobtools','create','-i',assembly_working,
                '-t',megablast_working,
                '-b', blobBAM,
                ]
    # could add -y spades not sure it matters though since providing -b bam
    if args.debug:
        logger.debug(blob_cmd)
    subprocess.run(blob_cmd,cwd=args.tmpdir)
    
#    blobtools view -i blobDB.json
#	grep -v '^#' blobDB.table.txt > blob.summary.txt
#	blobtools blobplot -i blobDB.json

	#calculate average coverage of largest 50 scaffolds
#	AVGCOV=$(egrep "$PHYLUM|no-hit" blob.summary.txt | head -n 50 | gawk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }')
#	MINCOV=$(bc <<< $AVGCOV*0.05)
#	MINCOV=${MINCOV%.*}
#	grep -v '^#' blob.summary.txt | gawk -v cov=$MINCOV '{ if($7 < cov) print $1}' > blob.low-coverage.txt

	#parse the blob summary output and filter for taxonomy and coverage, get rid of low coverage scaffolds, coverage less than 5% of average
#	grep -v '^#' blob.summary.txt | egrep -v "$PHYLUM|no-hit" | gawk '{print$1;}' > blob.contamination.txt

#        cat blob.contamination.txt blob.low-coverage.txt ncbi_remove-list.txt | sort | uniq > all.remove.txt

	#remove contaminant scaffolds
#	fasta_remove.py $SG1 all.remove.txt > spades.scafffolds.clean1.fasta
    
