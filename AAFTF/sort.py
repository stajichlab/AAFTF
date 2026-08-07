"""This module sorts FASTA sequences by size and renames headers."""

from Bio.SeqIO.FastaIO import SimpleFastaParser

from AAFTF.utility import softwrap, status


def run(parser, args):
    """Sort contig/scaffold file longest to shortest and rename."""
    status("Sorting sequences by length longest --> shortest")
    AllSeqs = {}
    with open(args.input) as fasta_in:
        for Header, Seq in SimpleFastaParser(fasta_in):
            if Header not in AllSeqs:
                if len(Seq) >= args.minlen:
                    AllSeqs[Header] = Seq
    sortSeqs = sorted(AllSeqs.items(), key=lambda item: len(item[1]), reverse=True)
    with open(args.out, "w") as fasta_out:
        for i, (Header, Seq) in enumerate(sortSeqs):
            fasta_out.write(f">{args.name}_{i + 1}\n{softwrap(Seq)}\n")

    status(f"Output written to: {args.out}")
