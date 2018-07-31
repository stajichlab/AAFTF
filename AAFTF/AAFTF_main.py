#!/usr/bin/env python3

# note structure of code taken from poretools https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import os.path
import sys
import argparse

import logging
logger = logging.getLogger('AAFTF')

# AAFTF imports
from AAFTF import version

def run_subtool(parser, args):
    if args.command == 'trim':
        import trim as submodule


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
    parser.add_argument("-v", "--version", help="Installed AAFTF version",
                        action="version",
                        version="%(prog)s " + str(AAFTF.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################

    ##########
    # combine
    ##########
    parser_trim = subparsers.add_parser('trim',
                                        help='Trim FASTQ input reads and Filter against known or defined contaminant library')
    # perhaps write this separately for singleton/unpaired read sets
    parser_trim.add_argument('forward',aliases=['fwd','1'],
                              dest='forward',
                              metavar='forward-reads',
                              required=True,
                              help='The name of the forward reads of paired-end FASTQ formatted reads.')
    parser_trim.add_argument('reverse',aliases=['rev','2'],
                            dest='reverse',
                            metavar='reverse-reads',
                            required=True,
                            help='The name of the reverse reads of paired-end FASTQ formatted reads.')

    parser_trim.set_defaults(func=run_subtool)

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
