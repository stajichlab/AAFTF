# trims fastq files of reads (typically Illumina)
# using trimmomatic or other specific trimmer (when written)
# attemps to remove vector and primer sequences

import sys, os, subprocess
from os.path import dirname
from AAFTF.utility import which_path
from AAFTF.utility import status
from AAFTF.utility import printCMD
from AAFTF.utility import Fzip_inplace
from AAFTF.utility import SafeRemove
from AAFTF.utility import getRAM
from AAFTF.utility import countfastq

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
    
    if not args.basename:
        if '_' in os.path.basename(args.left):
            args.basename = os.path.basename(args.left).split('_')[0]
        elif '.' in os.path.basename(args.left):
            args.basename = os.path.basename(args.left).split('.')[0]
        else:
            args.basename = os.path.basename(args.left)
    
    total = countfastq(args.left)
    if args.right:
        total = total*2
    status('Loading {:,} total reads'.format(total))
            
    DEVNULL = open(os.devnull, 'w')
    if args.method == 'bbduk':
        if args.memory:
            MEM='-Xmx{:}g'.format(args.memory)
        else:       
            MEM='-Xmx{:}g'.format(round(0.6*getRAM()))
            
        status('Adapter trimming using BBDuk')
        cmd = ['bbduk.sh', MEM, 'ref=adapters', 't={:}'.format(args.cpus), 'ktrim=r',
           'k=23', 'mink=11', 'minlen={:}'.format(args.minlength), 'hdist=1',
           'ftm=5', 'tpe', 'tbo', 'overwrite=true']
        if args.left and args.right:
            cmd += ['in1={:}'.format(args.left), 'in2={:}'.format(args.right), 
                    'out1={:}_1P.fastq.gz'.format(args.basename), 'out2={:}_2P.fastq.gz'.format(args.basename)]
        elif args.left:
            cmd += ['in={:}'.format(args.left), 'out={:}_1U.fastq.gz'.format(args.basename)]
        
        printCMD(cmd)
        if args.debug:
            subprocess.run(cmd)
        else:
            subprocess.run(cmd, stderr=DEVNULL)

        if args.right:
            clean = countfastq('{:}_1P.fastq.gz'.format(args.basename))
            clean = clean*2
            status('{:,} reads remaining and writing to file'.format(clean))
            status('Trimming finished:\n\tFor: {:}\n\tRev {:}'.format(
                        args.basename+'_1P.fastq.gz',
                        args.basename+'_2P.fastq.gz'))
            if not args.pipe:
                status('Your next command might be:\n\tAAFTF filter -l {:} -r {:} -o {:} -c {:}\n'.format(
                        args.basename+'_1P.fastq.gz', args.basename+'_2P.fastq.gz', args.basename, args.cpus))
        else:
            clean = countfastq('{:}_1U.fastq.gz'.format(args.basename))
            status('{:,} reads remaining and writing to file'.format(clean))
            status('Trimming finished:\n\tSingle: {:}'.format(
                        args.basename+'_1U.fastq.gz'))
            if not args.pipe:
                status('Your next command might be:\n\tAAFTF filter -l {:} -o {:} -c {:}\n'.format(
                        args.basename+'_1U.fastq.gz', args.basename, args.cpus))                 

    elif args.method == 'trimmomatic':
        #find path    
        trimmomatic_path = find_trimmomatic()
        if trimmomatic_path:
            jarfile = trimmomatic_path
        elif args.trimmomatic:
            jarfile = args.trimmomatic
        else:
            status('Trimmomatic cannot be found - please provide location of trimmomatic.jar file.')
            sys.exit(1)
        
        if jarfile:
            path_to_adaptors = args.trimmomatic_adaptors
            leadingwindow    = "LEADING:%d"%(args.trimmomatic_leadingwindow)
            trailingwindow   = "TRAILING:%d"%(args.trimmomatic_trailingwindow)
            slidingwindow    = "SLIDINGWINDOW:%s"%(args.trimmomatic_slidingwindow)

            quality = args.trimmomatic_quality
            quality = "-%s" % (quality) # add leading dash

            if not os.path.exists(path_to_adaptors):
                if args.right:
                    path_to_adaptors = dirname(jarfile)+"/adapters/TruSeq3-PE.fa"
                else:
                    path_to_adaptors = dirname(jarfile)+"/adapters/TruSeq3-SE.fa"
            
                if not os.path.exists(path_to_adaptors):
                    findpath=dirname(jarfile)
                    path_to_adaptors=""
                    while findpath:
                        if os.path.exists(findpath + "/share"):
                            if args.right:
                                path_to_adaptors=findpath+"/share/trimmomatic/adapters/TruSeq3-PE.fa"
                            else:
                                path_to_adaptors=findpath+"/share/trimmomatic/adapters/TruSeq3-SE.fa"
                            break
                        findpath=dirname(findpath)

                if not os.path.exists(path_to_adaptors):
                    status("Cannot find adaptors file, please specify manually")
                    status("Cannot find adaptors file, please specify manually")
                    return
        
            clipstr = args.trimmomatic_clip % (path_to_adaptors)
        
            cmd = []
        
            if args.left and args.right:
                cmd = ['java', '-jar', jarfile, 'PE',
                       '-threads',str(args.cpus),quality,
                       args.left,args.right,
                       args.basename+'_1P.fastq',
                       args.basename+'_1U.fastq',
                       args.basename+'_2P.fastq',
                       args.basename+'_2U.fastq',
                       clipstr, leadingwindow, trailingwindow,slidingwindow,
                       "MINLEN:%d" %(args.minlength) ]
            elif args.left and not args.right:
                cmd = ['java', '-jar', jarfile, 'SE',
                       '-threads',str(args.cpus),
                       quality,  args.left,
                       args.basename+'_1U.fastq',
                       clipstr, leadingwindow, trailingwindow,slidingwindow,
                       "MINLEN:%d" %(args.minlength) ]
            else:
                status("Must provide left and right pairs or single read set")
                return
        
            status('Running trimmomatic adapter and quality trimming')
            printCMD(cmd)
            if args.debug:
                subprocess.run(cmd)
            else:
                subprocess.run(cmd, stderr=DEVNULL)
            if args.right:
                status('Compressing trimmed PE FASTQ files')
                Fzip_inplace(args.basename+'_1P.fastq', args.cpus)
                Fzip_inplace(args.basename+'_2P.fastq', args.cpus)
                SafeRemove(args.basename+'_1U.fastq')
                SafeRemove(args.basename+'_2U.fastq')
                status('Trimming finished:\n\tFor: {:}\n\tRev {:}'.format(
                            args.basename+'_1P.fastq.gz',
                            args.basename+'_2P.fastq.gz'))
                if not args.pipe:
                    status('Your next command might be:\n\tAAFTF filter -l {:} -r {:} -o {:} -c {:}\n'.format(
                            args.basename+'_1P.fastq.gz', args.basename+'_2P.fastq.gz', args.basename, args.cpus)) 
            else:
                status('Compressing trimmed SE FASTQ file')
                Fzip_inplace(args.basename+'_1U.fastq', args.cpus)
                status('Trimming finished:\n\tSingle: {:}'.format(
                            args.basename+'_1U.fastq.gz'))
                if not args.pipe:
                    status('Your next command might be:\n\tAAFTF filter -l {:} -o {:} -c {:}\n'.format(
                            args.basename+'_1U.fastq.gz', args.basename, args.cpus))    

