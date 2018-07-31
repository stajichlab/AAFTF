#!/usr/bin/env python3

# note structure of code taken from poretools https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import os.path
import sys
import argparse

import logging
logger = logging.getLogger('AAFTF')

# AAFTF imports
from AAFTF import version
from AAFTF.version import __version__
myversion = __version__

def run_subtool(parser, args):
    if args.command == 'trim':
        import trim as submodule
    elif args.command == 'purge':
        import purge as submodule
    elif args.command == 'assemble':
        import assemble as submodule
    elif args.command == 'vecscreen':
        import vecscreen as submodule
    elif args.command == 'blobpurge':
        import blobpurge as submodule
    elif args.command == 'busco':
        import busco as submodule
    elif args.command == 'rmdup':
        import rmdup as submodule
    elif args.command == 'pilon':
        import pilon as submodule


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                            action="store_true",
                            dest="quiet")

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
    parser_trim = subparsers.add_parser('trim',
                                        help='Trim FASTQ input reads and Filter against known or defined contaminant library')
    # perhaps write this separately for singleton/unpaired read sets
    parser_trim.add_argument('--forward',type=str,nargs="+",
                              #aliases=['fwd','1'],
                              metavar='forward-reads',
                              required=True,
                              help='The name of the forward reads of paired-end FASTQ formatted reads.')

    parser_trim.add_argument('--reverse',type=str,nargs="+",
                            #aliases=['rev','2'],
                            metavar='reverse-reads',
                            required=True,
                            help='The name of the reverse reads of paired-end FASTQ formatted reads.')

    parser_trim.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                            help="Number of CPUs/threads to use.")elp="Temporary directory to store datafiles and processes in")
    ##########
    # purge
    ##########
    parser_trim.add_argument('-a','--screen_accession',type=str,nargs="+",
                            metavar='screening_accession',
                            required=False,
                            help="Genbank accession number(s) to screen out from initial reads.")
    parser_trim.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                            help="Number of CPUs/threads to use.")
    parser_trim.add_argument('--tmpdir',type=str,metavar='tmpdir',required=False,
                            help="Temporary directory to store datafiles and processes in")


    ### process args now ### 
    parser.set_defaults(func=run_subtool)
    args = parser.parse_args()
    
    if args.quiet:
        logger.setLevel(logging.ERROR)
    try:
        args.func(parser, args)

    except IOError as e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
