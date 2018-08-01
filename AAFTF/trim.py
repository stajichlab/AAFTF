import sys, os, subprocess

from os.path import dirname
#logging
import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):
    

    if args.trimmomatic:
        jarfile = args.trimmomatic
        path_to_adaptors=args.trimmomatic_adaptors
        clip    = args.trimmomatic_clip
        leadingwindow   = "LEADING:%d"%(args.trimmomatic_leadingwindow)
        trailingwindow  = "TRAILING:%d"%(args.trimmomatic_trailingwindow)
        slidingwindow   = "SLIDINGWINDOW:%s"%(args.trimmomatic_slidingwindow)
        minlen  = "MINLEN:%d" %(args.minlength)
        cpus    = args.cpus

        left    = args.left
        right   = args.right
        single  = args.single

        prefix  = args.prefix
        outdir  = args.outdir

        quality = args.trimmomatic_quality
        quality = "-%s" % (quality) # add leading dash

        if not path_to_adaptors:
            path_to_adaptors = dirname(jarfile)+"/adapters/TruSeq3-PE.fa"

        if not os.path.exists(path_to_adaptors):
            path_to_adaptors=dirname(dirname(jarfile))+"/share/trimmomatic/adapters/TruSeq3-PE.fa"

        if not os.path.exists(path_to_adaptors):
            print("Cannot find adaptors file, please specify manually")
            logger.info("Cannot find adaptors file, please specify manually")
            return
        
        clipstr = clip % (path_to_adaptors)

        
        cmd = []
        
        if left and right:
            cmd = ['java', '-jar', jarfile, 'PE',
                   '-threads',str(cpus),quality,
                   left,right,'-baseout', os.path.join(outdir,prefix),
                   clipstr, leadingwindow, trailingwindow,slidingwindow,minlen ]
        elif single:
            cmd = ['java', '-jar', jarfile, 'SE',
                   '-threads',str(cpus),
                   quality,  single,
                   '-baseout', os.path.join(outdir,prefix),
                   clipstr, leadingwindow, trailingwindow,slidingwindow,minlen ]
            

        logger.info("running cmd: %s" %(" ".join(cmd)))
        subprocess.call(cmd)

        # could change this up and support sickle or other processors
    else:
        print("Only trimmomatic supported as trim tool at the moment")
        logger.info("Only trimmomatic supported as trim tool at the moment")
