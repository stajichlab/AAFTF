# this module sorts FASTA sequences by size and renames headers
import sys
import os
import operator
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO

#logging
import logging
logger = logging.getLogger('AAFTF')

from AAFTF.utility import softwrap

def run(parser, args):
    logger.info('Sorting sequences by length longest --> shortest')
    AllSeqs = {}
    with open(args.input, 'rU') as fasta_in:
        for Header, Seq in SimpleFastaParser(fasta_in):
            if not Header in AllSeqs:
                AllSeqs[Header] = len(Seq)
    sortSeqs = sorted(AllSeqs.items(), key=operator.itemgetter(1),reverse=True)
    orderedSeqs = [i[0] for i in sortSeqs]
    SeqRecords = SeqIO.to_dict(SeqIO.parse(args.input, 'fasta'))
    with open(args.out, 'w') as fasta_out:
        for i,x in enumerate(orderedSeqs):
            fasta_out.write('>{:}_{:}\n{:}\n'.format(args.name, i+1, softwrap(str(SeqRecords[x].seq))))
    logger.info('Output written to: {:}'.format(args.out))
            



