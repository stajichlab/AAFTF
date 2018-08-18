# trims fastq files of reads (typically Illumina)
# using trimmomatic or other specific trimmer (when written)
# attemps to remove vector and primer sequences

import sys, os, subprocess

from os.path import dirname
#logging
import logging
logger = logging.getLogger('AAFTF')

from AAFTF.utility import which_path


# process trimming reads with trimmomatic
# Homebrew install of trimmomatic uses a shell script
'''
#!/bin/bash
exec java  -jar /usr/local/Cellar/trimmomatic/0.36/libexec/trimmomatic-0.36.jar "$@"
'''
#while bioconda install uses a python script that launches java apps

def find_trimmomatic():
    trim_path = which_path('trimmomatic')
    if trim_path:
        with open(os.path.abspath(which_path('trimmomatic')), 'rU') as trim_shell:
            firstLine = trim_shell.readline()
            if '#!/bin/bash' in firstLine: #then homebrew do routine to get jar location
                for line in trim_shell:         
                    if line.startswith('exec java'):
                        items = line.split(' ')
                        for x in items:
                            if x.endswith('.jar'):
                                return x
            elif '#!/usr/bin/env python' in firstLine:
                return os.path.join(os.path.dirname(os.path.realpath(trim_path)), 'trimmomatic.jar')
            else:
                return False         
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
    #find path    
    trimmomatic_path = find_trimmomatic()
    if trimmomatic_path:
        jarfile = trimmomatic_path
    elif args.trimmomatic:
        jarfile = args.trimmomatic
    else:
        logger.error('Trimmomatic cannot be found - please provide location of trimmomatic.jar file.')
        sys.exit(1)
        
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

