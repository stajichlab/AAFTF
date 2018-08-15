# run a set of default genome assembly using
# SPAdes. Additional tools could be supported or
# users may prefer to run their custom assembly and skip this step


import sys, os, subprocess

#logging
import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):
    
    if not args.outdir:
        args.outdir=args.indir

    spadescmd = ['spades.py','--threads',str(args.cpus),
                 '-k', '21,33,55,77,99,127','--mem',args.memory,'--careful',
                 '-o',os.path.join(args.outdir,args.prefix) ]
    if args.pairing:
        left = os.path.join(args.indir,args.prefix + "_cleaned_1.fq.gz")
        right = os.path.join(args.indir,args.prefix + "_cleaned_2.fq.gz")
        
        if not os.path.exists(left):
            interleaved = os.path.join(args.indir,args.prefix + "_cleaned_12.fq.gz")
            if os.path.exists(interleaved):
                spadescmd.append('-12')
                spadescmd.append(interleaved)                
            else:
                logger.info("no interleaved file %s and no left file %s" %(interleaved,left))
                print("no interleaved file %s and no left file %s" %(interleaved,left))
                return
        else:
            spadescmd.append('--pe1-1')
            spadescmd.append(left)
            spadescmd.append('--pe1-2')
            spadescmd.append(right)

    else:
        left = os.path.join(args.indir,args.prefix + "_cleaned.1.gz")
        if not os.path.exists(left):
            print("No left/single file found for unpaired run of asm: %s" % (left))
            return
        
        spadescmd.append('-s')
        spadescmd.append(left)

    # now run the spades job
    print(spadescmd)
    subprocess.run(spadescmd)
