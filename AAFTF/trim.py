import sys, os, subprocess

from os.path import dirname
#logging
import logging
logger = logging.getLogger('AAFTF')

from AAFTF.utility import which_path


# process trimming reads with trimmomatic

def find_trimmomatic():
    if which_path('trimmomatic'):
        with open(os.path.abspath(which_path('trimmomatic')), 'rU') as trim_shell:
            for line in trim_shell:
                if line.startswith('exec java'):
                    items = line.split(' ')
                    for x in items:
                        if x.endswith('.jar'):
                            print(x)
                            return x
                        
    else:
        return False
    
def run(parser,args):    

    if not args.outdir:
        args.outdir = dirname(args.left)
    os.makedirs(args.outdir,exist_ok=True)

    left_expected = os.path.join(args.outdir,args.prefix)+"_1P.fastq"
    if ( os.path.exists(left_expected) and
         os.path.getctime(left_expected) >
         os.path.getctime(args.left) ):
        logger.info("Already ran trimming on %s %s" % (args.left,args.right))
        return
    
    if args.trimmomatic:
        jarfile      = args.trimmomatic
    else:
        jarfile      = find_trimmomatic()
        
    if jarfile:
        path_to_adaptors = args.trimmomatic_adaptors
        leadingwindow    = "LEADING:%d"%(args.trimmomatic_leadingwindow)
        trailingwindow   = "TRAILING:%d"%(args.trimmomatic_trailingwindow)
        slidingwindow    = "SLIDINGWINDOW:%s"%(args.trimmomatic_slidingwindow)

        quality = args.trimmomatic_quality
        quality = "-%s" % (quality) # add leading dash

        if not os.path.exists(path_to_adaptors):
            path_to_adaptors = dirname(jarfile)+"/adapters/TruSeq3-PE.fa"
            
            if not os.path.exists(path_to_adaptors):
                findpath=dirname(jarfile)
                path_to_adaptors=""
                while findpath:
                    if os.path.exists(findpath + "/share"):
                        path_to_adaptors=findpath+"/share/trimmomatic/adapters/TruSeq3-PE.fa"
                        break
                    findpath=dirname(findpath)

            if not os.path.exists(path_to_adaptors):
                print("Cannot find adaptors file, please specify manually")
                logger.info("Cannot find adaptors file, please specify manually")
                return
        
        clipstr = args.trimmomatic_clip % (path_to_adaptors)

        
        cmd = []
        
        if args.left and args.right:
            cmd = ['java', '-jar', jarfile, 'PE',
                   '-threads',str(args.cpus),quality,
                   args.left,args.right,
                   os.path.join(args.outdir,args.prefix+'_1P.fastq'),
                   os.path.join(args.outdir,args.prefix+'_1U.fastq'),
                   os.path.join(args.outdir,args.prefix+'_2P.fastq'),
                   os.path.join(args.outdir,args.prefix+'_2U.fastq'),
                   clipstr, leadingwindow, trailingwindow,slidingwindow,
                   "MINLEN:%d" %(args.minlength) ]
        elif args.single:
            cmd = ['java', '-jar', jarfile, 'SE',
                   '-threads',str(args.cpus),
                   quality,  args.single,
                   os.path.join(args.outdir,args.prefix+'_1U.fastq'),
                   clipstr, leadingwindow, trailingwindow,slidingwindow,
                   "MINLEN:%d" %(args.minlength) ]
        else:
            logger.error("Must provide left and right pairs or single read set")
            return
        
        logger.info("running cmd: %s" %(" ".join(cmd)))
        subprocess.run(cmd)

    else:
        print("Only trimmomatic or sickle supported as trim tool at the moment")
        logger.info("Only trimmomatic or sickle supported as trim tool at the moment")
