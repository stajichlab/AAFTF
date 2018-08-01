#!/usr/bin/env python3

# note structure of code taken from poretools https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import os.path
import sys
import argparse

import logging
logger = logging.getLogger('AAFTF')

# AAFTF imports
from AAFTF.version import __version__
myversion = __version__

def run_subtool(parser, args):
    if args.command == 'trim':
        import AAFTF.trim as submodule
    elif args.command == 'filter':
        import AAFTF.filter as submodule
    elif args.command == 'assemble':
        import AAFTF.assemble as submodule
    elif args.command == 'vecscreen':
        import AAFTF.vecscreen as submodule
    elif args.command == 'blobfilter':
        import AAFTF.blobfilter as submodule
    elif args.command == 'busco':
        import AAFTF.busco as submodule
    elif args.command == 'rmdup':
        import AAFTF.rmdup as submodule
    elif args.command == 'pilon':
        import AAFTF.pilon as submodule
    else:
        parser.parse_args('')
        return
    # run the chosen submodule.
    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                            action="store_true",
                            dest="quiet")
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    logging.basicConfig()

    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='AAFTF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")
    parser.add_argument("-v", "--version", help="Installed AAFTF version",
                        action="version",
                        version="%(prog)s " + str(myversion))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)
    #########################################
    # create the individual tool parsers
    #########################################

    ##########
    # trim
    ##########
    # arguments
    # --trimmomatic or --sickle: arguments are path to JAR or application respectively
    # assume java is PATH already for trimmomatic
    # -c / --cpus: CPUs [default=1]
    # -o / --outdir: write outdir
    # -p / --prefix: outfile prefix
    # -ml / --minlength: min read length
    
    # --tmpdir: temporary directory to write tmp files [optional]
    # read info, either paired data are required or singleton
    # --left: left or forward reads
    # --right: right or reverse reads
    # currently singleton / unpaired reads not supported?
    
    parser_trim = subparsers.add_parser('trim',
       description="This comamnd trims reads in FASTQ format to remove low quality reads and trim adaptor sequences",
       help='Trim FASTQ input reads')
    
    parser_trim.add_argument('-p','--prefix',type=str,
                             required=True,
                             help="Output Prefix")

    parser_trim.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                              help="Number of CPUs/threads to use.")
    
    parser_trim.add_argument('--tmpdir',type=str,
                             required=False,
    help="Temporary directory to store datafiles and processes in")

    parser_trim.add_argument('-o','--outdir',type=str,
                             required=False,
    help="Output directory for trimmed reads")

    parser_trim.add_argument('-ml','--minlength',type=int,
                             default=75,
                             required=False,
                             help="Minimum read length after trimming, default: 75")

    tool_group = parser_trim.add_mutually_exclusive_group(required=True)

    tool_group.add_argument('--trimmomatic','--jar', metavar='trimmomatic_jar',
                            type=str,required=False,
                            help='Trimmomatic JAR path')
    
    tool_group.add_argument('--sickle',type=str,
                            required=False,
                            nargs="?",
                            const='1',
                            default='0',
            help="Use sickle program for read processing (specify path to prog if not in PATH already)")
    
    trimmomatic_group = parser_trim.add_argument_group(title='Trimmomatic options',
                                              description="Trimmomatic trimming options")

    trimmomatic_group.add_argument('--trimmomatic_adaptors',
                                   default="TruSeq3-PE.fa",
                                   help="Trimmomatic adaptor file, default: TruSeq3-PE.fa")

    trimmomatic_group.add_argument('--trimmomatic_clip',
                                   default="ILLUMINACLIP:%s:2:30:10",
                                   help="Trimmomatic clipping, default: ILLUMINACLIP:TruSeq3-PE.fa:2:30:10")

    trimmomatic_group.add_argument('--trimmomatic_leadingwindow',
                                   default="3",type=int,
                                   help="Trimmomatic window processing arguments, default: LEADING:3")

    trimmomatic_group.add_argument('--trimmomatic_trailingwindow',
                                   default="3",type=int,
                                   help="Trimmomatic window processing arguments, default: TRAILING:3")

    trimmomatic_group.add_argument('--trimmomatic_slidingwindow',
                                   default="4:15",type=str,
                                   help="Trimmomatic window processing arguments, default: SLIDINGWINDOW:4:15")

    
    paired_reads = parser_trim.add_argument_group(title='Paired Reads',
                                                  description="Paired Read FASTQ files")
    
    trimmomatic_group.add_argument('--trimmomatic_quality',
                                   default="phred33",
                                   help="Trimmomatic quality encoding -phred33 or phred64")

    paired_reads.add_argument('--left',type=str,
                              required=False,
            help='The name of the left/forward reads of paired-end FASTQ formatted reads.')

    paired_reads.add_argument('--right',type=str,
                              required=False,
            help='The name of the right/reverse reads of paired-end FASTQ formatted reads.')

    # perhaps write this separately for singleton/unpaired read sets
    parser_trim.add_argument('--single',type=str,
                             required=False,
    help='The name of the reverse reads of paired-end FASTQ formatted reads.')

    ##########
    # filter
    ##########
    # arguments
    # -i / --indir:  input dir
    # -o / --outdir: write outdir
    # -p / --prefix: outfile prefix

    parser_filter = subparsers.add_parser('filter',
       description="Filter reads which match contaminant databases such as phiX",
       help='Filter contaminanting reads')

    parser_filter.add_argument('-a','--screen_accessions',type=str,
                               metavar='accessions',nargs="*",
                            help="Genbank accession number(s) to screen out from initial reads.")
    parser_filter.add_argument('-c','--cpus',type=int,metavar="cpus",default=1,
                              help="Number of CPUs/threads to use.")
    parser_filter.add_argument('--tmpdir',type=str,metavar='tmpdir',
                              help="Temporary directory to store datafiles and processes in")

    parser_filter.add_argument('-i','--indir',type=str,
                             required=True,
    help="Directory for input of trimmed reads")

    parser_filter.add_argument('-o','--outdir',type=str,
                             required=False,
                             help="Directory for filtered reads (defaults to indir)")

    tool_group = parser_filter.add_mutually_exclusive_group(required=True)

    tool_group.add_argument('--bowtie2',
                            type=str,required=False,default='1',const='1',nargs='?',
                            help='Bowtie2 executable path (specify path to bowtie2 if not in PATH already)')
    
    tool_group.add_argument('--bwa',type=str,
                            required=False,
                            nargs="?",
                            const='1',
                            default='0',
            help="Use bwa for read filtering (specify path to bwa prog if not in PATH already)")

    tool_group.add_argument('--bbmap',type=str,
                            required=False,
                            nargs="?",
                            const='1',
                            default='0',
            help="Use bbmap for read filtering (specify path to bbmap prog if not in PATH already)")


    
    parser.set_defaults(func=run_subtool)

    ### process args now ###
    # if no args then print help and exit
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()
    fh = logging.FileHandler('AAFTF.log')
    fh.setLevel(logging.INFO)
    
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    
    if args.quiet:
        logger.setLevel(logging.ERROR)
    try:
        args.func(parser, args)

    except IOError as e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
