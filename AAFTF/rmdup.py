import sys, os

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

    input = args.in # this is required
    out = args.out     # this is also required

    print(sys.argv[0])
    cmd = ['scripts/funannotate-contig_cleaner.py',
           '--input',input,
           '--out',out]

    if args.percent_id:
        cmd.append('--pident')
        cmd.append(args.percent_id)

    print(cmd)
           
           

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)

    if args.percent_id:
        percentid_cutoff = args.percent_id
    else:
        percentid_cutoff = default_percent_id_cutoff
        
