import sys, os

#logging
import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):

    if parser.tool == "trimmomatic":
        jarfile = "trimmomatic.jar" # where to  get this globally set
        path_to_adaptors=os.path.dirname(jarfile)+"/adapters/TruSeq3-PE.fa"

        #outfwd,outrev,outfwd_unpair,outrev_unpair
        cmds=['ILLUMINACLIP:%s:2:30:10'%(path_to_adaptors),'LEADING:3','TRAILING:3','SLIDINGWINDOW:4:15','MINLEN:75']
        ['java','-jar', jarfile,'PE', args.forward,args.reverse,
            outfwd,outrev,outfwd_unpair,outrev_unpair,
            cmds ]
        logger.info("Ran Trimmomatic with files %s %s.\n" % \
                (args.forward, args.reverse))

    else:
        # could change this up and support sickle or other processors
        print("Only trimmomatic supported as trim tool at the moment")
        logger.info("Only trimmomatic supported as trim tool at the moment")
