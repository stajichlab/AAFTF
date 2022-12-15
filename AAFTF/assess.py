"""Assess quality statistics of a genome assembly.

This simply gives GC%, N50, L50, Min, Max
contig statistics.
"""
import gzip
import os
import re

from Bio import SeqIO

from AAFTF.utility import status


def genome_asm_stats(fasta_file, output_handle, telomere_repeat, n_minimum):
    """Calculate genome assembly statistics."""
    lengths = []
    # could be smart here and handle compressed files?
    GC = 0
    if fasta_file.endswith(".gz"):
        seqio = SeqIO.parse(gzip.open(fasta_file, mode="rt"), "fasta")
    else:
        seqio = SeqIO.parse(fasta_file, "fasta")

    telomere_stats = {'FORWARD TELOMERES': 0,
                      'REVERSE TELOMERES': 0,
                      'BOTH TELOMERES': 0}
    for record in seqio:
        lengths.append(len(record))
        forward, reverse = findTelomere(record.seq, telomere_repeat, n_minimum)
        if forward:
            telomere_stats['FORWARD TELOMERES'] += 1
        if reverse:
            telomere_stats['REVERSE TELOMERES'] += 1
        if forward and reverse:
            telomere_stats['BOTH TELOMERES'] += 1

        GC += sum(record.seq.count(x) for x in ["G", "C", "g", "c", "S", "s"])

    lengths.sort()
    total_len = sum(lengths)
    GC = 100.0 * (GC / total_len)
    l50 = 0
    n50 = 0
    l90 = 0
    n90 = 0
    cumulsum = 0
    i = 1
    for n in reversed(lengths):
        cumulsum += n
        if n50 == 0 and cumulsum >= total_len * 0.5:
            n50 = n
            l50 = i
        if n90 == 0 and cumulsum >= total_len * 0.9:
            n90 = n
            l90 = i

        i += 1
    report = "Assembly statistics for: %s\n" % (fasta_file)
    report += "%15s  =  %d\n" % ('CONTIG COUNT', len(lengths))
    report += "%15s  =  %d\n" % ('TOTAL LENGTH', total_len)
    report += "%15s  =  %d\n" % ('MIN', lengths[0])
    report += "%15s  =  %d\n" % ('MAX', lengths[-1])
    report += "%15s  =  %d\n" % ('MEDIAN', lengths[int(len(lengths)/2)])
    report += "%15s  =  %.2f\n" % ('MEAN',  total_len/len(lengths))
    report += "%15s  =  %d\n" % ('L50',  l50)
    report += "%15s  =  %d\n" % ('N50',  n50)
    report += "%15s  =  %d\n" % ('L90',  l90)
    report += "%15s  =  %d\n" % ('N90',  n90)
    report += "%15s  =  %.2f\n" % ('GC%', GC)
    for f in sorted(telomere_stats):
        report += "%15s  =  %.2f\n" % (f, telomere_stats[f])

    print(report)
    if output_handle:
        output_handle.write(report)


def revcomp(seq):
    """Reverse complement sequence or regexp.

    Using code based on find_telomeres.py from Markus Hiltunen.
    https://github.com/markhilt/genome_analysis_tools
    """
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
        else:  # At this point we don't care about IUPAC coded bases
            revcomped_seq.append("N")
    return "".join(revcomped_seq)


def findTelomere(seq, monomer, n):
    """Takes nucleotide sequence and checks if the sequence contains telomere repeats.

    Using code based on find_telomeres.py from Markus Hiltunen.
    https://github.com/markhilt/genome_analysis_tools
    """
    # Look within first and last 200 bp for repeats
    start = seq[:100].upper()
    end = seq[-100:].upper()

    forward, reverse = False, False

    # Look for the monomer repeat n number of times.
    if re.search(monomer*n, start):
        forward = True
    rev_monomer = revcomp(monomer)
    if re.search(rev_monomer*n, end):
        reverse = True

    return forward, reverse


def run(parser, args):
    """This is the general run command to calculate the genome statistics.

    This function will also attempt to find the telomere repeats and count these.
    """
    if not os.path.exists(args.input):
        status("Inputfile %s was not readable, check parameters" % (
            args.input))

    output_handle = None

    if args.report:
        output_handle = open(args.report, "w")

    genome_asm_stats(args.input, output_handle, args.telomere_monomer, args.telomere_n_repeat)
