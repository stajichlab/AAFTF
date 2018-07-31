import sys, os

#import os.path.dirname as dirname
#logging
import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):
    

    if args.trimmomatic:
        jarfile = args.trimmomatic
        path_to_adaptors=args.trimmomatic_adaptors

        if not path_to_adaptors:
            path_to_adaptors = os.path.dirname(jarfile)+"/adapters/TruSeq3-PE.fa"
        print(jarfile)
        if not os.path.exists(path_to_adaptors):
            path_to_adaptors=os.path.dirname(jarfile)+"../share/trimmomatic/adapters/TruSeq3-PE.fa"
            
        if not os.path.exists(path_to_adaptors):
            logger.info("Cannot find adaptors file, please specify manually")
            
        #        logger.info("Ran Trimmomatic with files %s %s.\n" % \
            #                (args.forward, args.reverse))

        # could change this up and support sickle or other processors
    else:
        print("Only trimmomatic supported as trim tool at the moment")
        logger.info("Only trimmomatic supported as trim tool at the moment")
