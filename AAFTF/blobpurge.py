import sys, os
import shutil

from subprocess import call, Popen, PIPE, STDOUT

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#logging
import logging
logger = logging.getLogger('AAFTF')


def run(parser,args):

    if not args.tmpdir:
        args.tmpdir = 'working_AAFTF'

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)

    # Phylum
