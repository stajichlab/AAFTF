#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Reads a fasta file and uses regex to find telomere repeats.

Copyright (c) 2021, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
'''

import argparse
import re

__version__ = "0.2"

parser = argparse.ArgumentParser(description="Reads a fasta file and uses regex to find telomere repeats.")
parser.add_argument("reference", help="Reference genome fasta file. Required.", type = str)
parser.add_argument("-m", "--monomer", \
                    help="Custom monomer to look for. Run separately for reverse telomeres. [TAA[C]+]", \
                    type = str, default = "TAA[C]+")
parser.add_argument("-n", "--n_repeats", \
                    help="Minimum number of monomer repeats. [2]", \
                    type = int, default = 2)
parser.add_argument("-v","--version", help="Print version and quit.", \
                    action = "version", \
                    version = "find_telomeres v.{}".format(__version__))
args = parser.parse_args()

def readFasta():
    fa = {}
    with open(args.reference, "r") as fasta:
        sequence = None
        for line in fasta:
            line = line.strip()

            if line.startswith(">") and sequence == None:
                header = line[1:]
                sequence = []

            elif line.startswith(">") and sequence != None:
                # If new fasta entry, add old one to dict, pick new header and reset sequence
                fa[header] = "".join(sequence)
                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

        # Last passthrough won't have any new entries, just add the remaining sequence
        fa[header] = "".join(sequence)

    return fa

def revcomp(seq):
    ''' Reverse complement sequence. Also the regex. Kinda.
    '''
    revcomped_seq = []
    for nucl in seq[::-1]:
        if nucl == "A" or nucl == "a":
            revcomped_seq.append("T")
        elif nucl == "T" or nucl == "t":
            revcomped_seq.append("A")
        elif nucl == "C" or nucl == "c":
            revcomped_seq.append("G")
        elif nucl == "G" or nucl == "g":
            revcomped_seq.append("C")
        elif nucl == "[":
            revcomped_seq.append("]+")
        elif nucl == "]":
            revcomped_seq.append("[")
        elif nucl == "+":
            continue
        else: # At this point we don't care about IUPAC coded bases
            revcomped_seq.append("N")
    return "".join(revcomped_seq)

def findTelomere(seq, monomer, n):
    '''
    Takes nucleotide sequence and checks if the sequence contains telomere repeats.
    '''
    # Look within first and last 200 bp for repeats
    start = seq[:100].upper()
    end = seq[-100:].upper()

    forward, reverse = False, False
    '''
    # Look for TAACCCC... and TTAGGGG...
    if re.search("TAA[C]+TAA[C]+", start):
        forward = True
    if re.search("TTA[G]+TTA[G]+", end):
        reverse = True
    '''

    # Look for the monomer repeat n number of times.
    if re.search(monomer*n, start):
        forward = True
    rev_monomer = revcomp(monomer)
    if re.search(rev_monomer*n, end):
        reverse = True

    return forward, reverse

def main():
    fasta = readFasta()
    n_f, n_r = 0,0

    for header, seq in fasta.items():
        forward, reverse = findTelomere(seq,args.monomer,args.n_repeats)

        if forward == True:
            print("{}\tforward\t{}".format(header, seq[:100]))
            n_f += 1
        if reverse == True:
            print("{}\treverse\t{}".format(header, seq[-100:]))
            n_r += 1

    print("\nTelomeres found: {} ({} forward, {} reverse)".format(str(n_f+n_r),n_f,n_r))

if __name__ == "__main__":
    main()
