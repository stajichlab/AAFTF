import sys, os, subprocess

from os.path import dirname
#logging
import logging
logger = logging.getLogger('AAFTF')

# process trimming reads with trimmomatic or sickle

def run(parser,args):    

    if not args.outdir:
        args.outdir = dirname(args.left)
    os.makedirs(args.outdir,exist_ok=True)

    left_expected = os.path.join(args.outdir,args.prefix)+"_1P"
    if ( os.path.exists(left_expected) and
         os.path.getctime(left_expected) >
         os.path.getctime(args.left) ):
        logger.info("Already ran trimming on %s %s" % (args.left,args.right))
        return
            
    if args.trimmomatic:
        jarfile          = args.trimmomatic
        path_to_adaptors = args.trimmomatic_adaptors
        leadingwindow    = "LEADING:%d"%(args.trimmomatic_leadingwindow)
        trailingwindow   = "TRAILING:%d"%(args.trimmomatic_trailingwindow)
        slidingwindow    = "SLIDINGWINDOW:%s"%(args.trimmomatic_slidingwindow)

        quality = args.trimmomatic_quality
        quality = "-%s" % (quality) # add leading dash

        if not path_to_adaptors:
            path_to_adaptors = os.path.join(dirname(jarfile),
                                            "/adapters/TruSeq3-PE.fa")
            
        logger.info("Tried to find adaptors in %s"%(path_to_adaptors))

        if not os.path.exists(path_to_adaptors):
            path_to_adaptors=os.path.join(dirname(dirname(jarfile)),
                                          "/share/trimmomatic/adapters/TruSeq3-PE.fa")

        logger.info("Tried to find adaptors in %s"%(path_to_adaptors))
        if not os.path.exists(path_to_adaptors):
            print("Cannot find adaptors file, please specify manually")
            logger.info("Cannot find adaptors file in %s, please specify manually" % (path_to_adaptors))
            return
        
        clipstr = args.trimmomatic_clip % (path_to_adaptors)

        
        cmd = []
        
        if args.left and args.right:
            cmd = ['java', '-jar', jarfile, 'PE',
                   '-threads',str(args.cpus),quality,
                   args.left,args.right,'-baseout',
                   os.path.join(args.outdir,args.prefix),
                   clipstr, leadingwindow, trailingwindow,slidingwindow,
                   "MINLEN:%d" %(args.minlength) ]
        elif args.single:
            cmd = ['java', '-jar', jarfile, 'SE',
                   '-threads',str(args.cpus),
                   quality,  args.single,
                   '-baseout', os.path.join(args.outdir,args.prefix),
                   clipstr, leadingwindow, trailingwindow,slidingwindow,
                   "MINLEN:%d" %(args.minlength) ]
        else:
            logger.error("Must provide left and right pairs or single read set")
            return

        logger.info("running cmd: %s" %(" ".join(cmd)))
        subprocess.run(cmd)

        # could change this up and support sickle or other processors
    elif args.sickle:
        cmd = []

        if args.left and args.right:
            cmd = ['sickle', 'pe', '-t', 'sanger',
                   '-f',args.left,'-r',args.right,'-l',str(args.minlength),
                   '-o',os.path.join(args.outdir,args.prefix)+"_1P",
                   '-p',os.path.join(args.outdir,args.prefix)+"_2P",
                   '-s',os.path.join(args.outdir,args.prefix)+"_1U"
            ]
        elif args.single:
            cmd = ['sickle', 'se', '-t', 'sanger',
                   '-f',single, '-l',args.minlength,
                   '-o',os.path.join(args.outdir,args.prefix)+"_1U" ]
        else:
            logger.error("Must provide left and right pairs or single read set")
            return
        print(cmd)
        logger.info("running cmd: %s" %(" ".join(cmd)))
        subprocess.run(cmd)


    else:
        print("Only trimmomatic or sickle supported as trim tool at the moment")
        logger.info("Only trimmomatic or sickle supported as trim tool at the moment")
