import sys, os
import urllib.request

#logging
import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):

    
    outdir = args.outdir

    if args.bowtie2:

    elif args.bwa:

    elif args.bbmap:

    else:
        print("Must specify bowtie2,bwa, or bbmap for filtering")
        logger.info("Must specify bowtie2,bwa, or bbmap for filtering")

            
