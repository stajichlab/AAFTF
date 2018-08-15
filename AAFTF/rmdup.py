import sys, os

from subprocess import call, Popen, PIPE, STDOUT

#logging
import logging
logger = logging.getLogger('AAFTF')

def run(parser,args):

    # here we will just run funannotate clean
    # note that this generally runs in the current working dir and I
    # thing we need to fix it to better handle running multiple analysis
    # from same working dir - so a prefix or tempfolder is needed
    # maybe an update to that code

    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)

    input = args.input
    out = args.out     # this is also required

    print(sys.argv[0])
    cmd = ['scripts/funannotate-contig_cleaner.py',
           '--input',input,
           '--out',out]

    if args.percent_id:
        cmd.append('--pident')
        cmd.append(str(args.percent_id))

    if args.coverage:
        cmd.append('--coverage')
        cmd.append(str(args.coverage))
        
    if args.minlen:
        cmd.append('--minlen')
        cmd.append(str(args.minlen))

    if args.exhaustive:
        cmd.append('--exhaustive')

    if args.method:
        cmd.append('--method')
        cmd.append(args.method)

    if args.debug:
        cmd.append('--debug')

    call(cmd)
           

        
