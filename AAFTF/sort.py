"""This module sorts FASTA sequences by size and renames headers."""
import operator

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from AAFTF.utility import softwrap, status


def run(parser, args):
    """Sort contig/scaffold file longest to shortest and rename."""
    status('Sorting sequences by length longest --> shortest')
    AllSeqs = {}
    with open(args.input) as fasta_in:
        for Header, Seq in SimpleFastaParser(fasta_in):
            if Header not in AllSeqs:
                if len(Seq) >= args.minlen:
                    AllSeqs[Header] = len(Seq)
    sortSeqs = sorted(
        AllSeqs.items(), key=operator.itemgetter(1), reverse=True)
    orderedSeqs = [i[0] for i in sortSeqs]
    SeqRecords = SeqIO.to_dict(SeqIO.parse(args.input, 'fasta'))
    with open(args.out, 'w') as fasta_out:
        for i, x in enumerate(orderedSeqs):
            fasta_out.write('>{:}_{:}\n{:}\n'.
                            format(args.name, i+1,
                                   softwrap(str(SeqRecords[x].seq))))

    status(f'Output written to: {args.out}')
