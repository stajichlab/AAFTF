# assessment of assembly - this can be as simple as
# N50/L50/MAX contig stats

import sys, os
from Bio import SeqIO
from AAFTF.utility import status


def genome_asm_stats(fasta_file,output_handle):

    lengths = []
    # could be smart here and handle compressed files?
    GC  = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record))
        GC += sum(record.seq.count(x) for x in ["G", "C", "g", "c", "S", "s"])

    lengths.sort()
    total_len = sum(lengths)
    GC  = 100.0 * (GC / total_len)
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
    report += "%15s  =  %d\n" % ('CONTIG COUNT',len(lengths))
    report += "%15s  =  %d\n" % ('TOTAL LENGTH',total_len)
    report += "%15s  =  %d\n" % ('MIN',lengths[0])
    report += "%15s  =  %d\n" % ('MAX',lengths[-1])
    report += "%15s  =  %d\n" % ('MEDIAN',lengths[int(len(lengths)/2)])
    report += "%15s  =  %.2f\n" % ('MEAN',  total_len/len(lengths))
    report += "%15s  =  %d\n" % ('L50',  l50)
    report += "%15s  =  %d\n" % ('N50',  n50)
    report += "%15s  =  %d\n" % ('L90',  l90)
    report += "%15s  =  %d\n" % ('N90',  n90)
    report += "%15s  =  %.2f\n" % ('GC%', GC)

    print(report)
    if output_handle:
        output_handle.write(report)

def run(parser,args):

    if not os.path.exists(args.input):
        status("Inputfile %s was not readable, check parameters" % (args.input))

    output_handle=None

    if args.report:
        output_handle = open(args.report,"w")

    genome_asm_stats(args.input, output_handle)
