#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
import uuid
import shutil
from Bio.SeqIO.FastaIO import SimpleFastaParser
from AAFTF.utility import softwrap, execute, RevComp, GuessRL, getRAM, printCMD, status, which_path


def orient_to_start(fasta_in, fasta_out, folder='.', start=False):
    # if not starting, then use cytochrome oxidase (cob)
    startFile = os.path.join(folder, '{}.fasta'.format(uuid.uuid4()))
    if not start:
        # generated as spoa consensus from select fungal cob genes
        # move this to a configurable file
        cob1 = 'atgagaattttaaaaagtcatcctttattaaaattagttaatagttatattattgattcaccacaaccttctaatattagttatttatgaaattttggatctttattagctttatgtttagttatacaaattgtaactggtgttacattagctatgcactatacacctaatgttgatttagcttttaattctgtagaacatattatgagagatgtaaataatggttgattaataagatatttacatgctaatactgcttcagcattctttttcttagttatatttacatataggtagaggattatattatggttcatataaatcacctagaacattaacatgagctattgg'
        with open(startFile, 'w') as outfile:
            outfile.write('>COB\n{}\n'.format(softwrap(cob1)))
    else:
        shutil.copyfile(start, startFile)

    # load sequence into dictionary
    initial_seq = ''
    header = ''
    with open(fasta_in, 'r') as infile:
        for title, seq in SimpleFastaParser(infile):
            initial_seq = seq
            header = title

    alignments = []
    minimap2_cmd = ['minimap2', '-x', 'map-ont', '-c', fasta_in, startFile]
    for line in execute(minimap2_cmd, '.'):
        cols = line.rstrip().split('\t')
        alignments.append(cols)
    if len(alignments) == 1:
        ref_strand = cols[4]
        ref_offset = int(cols[2])
        if ref_strand == '-':
            ref_start = int(cols[8]) + ref_offset
        else:
            ref_start = int(cols[7]) - ref_offset
        rotated = initial_seq[ref_start:] + initial_seq[:ref_start]
        if ref_strand == '-':
            rotated = RevComp(rotated)
        with open(fasta_out, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format('mt', softwrap(rotated)))
    elif len(alignments) == 0:
        status('ERROR: unable to rotate because did not find --starting sequence\n')
        with open(fasta_out, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format('mt', softwrap(initial_seq)))
    elif len(alignments) > 1:
        status('ERROR: unable to rotate because found multiple alignments\n')
        for x in alignments:
            sys.stderr.write('{}\n'.format(x))
        with open(fasta_out, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format('mt', softwrap(initial_seq)))
    if os.path.isfile(startFile):
        os.remove(startFile)

def run(parser,args):
    # first check if NOVOplasty and minimap2 are installed, else exit
    programs = ['NOVOplasty.pl', 'minimap2']
    for x in programs:
        if not which_path(x):
            status('ERROR: {} is not installed, exiting'.format(x))
            sys.exit(1)
    # first we need to generate working directory
    unique_id = str(uuid.uuid4())[:8]
    if not args.workdir:
        args.workdir = 'mito_'+unique_id
    if not os.path.isdir(args.workdir):
        os.makedirs(args.workdir)

    # now estimate read lengths of FASTQ
    read_len = GuessRL(args.left)

    # check for seed sequence, otherwise write one
    if not args.seed:
        if not args.reference:
            seedFasta = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data','mito-seed.fasta'))
        else:
            seedFasta = os.path.abspath(args.reference)
    else:
        seedFasta = os.path.abspath(args.seed)

    # now write the novoplasty config file
    defaultConfig = os.path.join(os.path.dirname(__file__), 'data','novoplasty-config.txt')
    novoConfig = os.path.join(args.workdir, 'novo-config.txt')
    if args.reference:
        refgenome = os.path.abspath(args.reference)
    else:
        refgenome = ''
    checkWords = ("<PROJECT>", "<MINLEN>", "<MAXLEN>", "<MAXMEM>",
                  "<SEED>", "<READLEN>", "<FORWARD>", "<REVERSE>",
                  "<REFERENCE>")
    repWords = (unique_id, str(args.minlen), str(args.maxlen), str(int(getRAM()*.75)),
                seedFasta,  str(read_len), os.path.abspath(args.left), os.path.abspath(args.right),
                refgenome)
    with open(novoConfig, 'w') as outfile:
        with open(defaultConfig, 'r') as infile:
            for line in infile:
                for check, rep in zip(checkWords, repWords):
                    line = line.replace(check, rep)
                outfile.write(line)

    # now we can finally run NOVOplasty.pl
    status('De novo assembling mitochondrial genome using NOVOplasty')
    cmd = ['NOVOPlasty.pl', '-c', 'novo-config.txt']
    printCMD(cmd)
    novolog = os.path.join(args.workdir, 'novoplasty.log')
    with open(novolog, 'w') as logfile:
        p1 = subprocess.Popen(cmd, cwd=args.workdir, stdout=logfile, stderr=logfile)
        p1.communicate()

    # now parse the results
    draftMito = None
    circular = False
    for f in os.listdir(args.workdir):
        if f.startswith('Circularized_assembly_'):
            draftMito = os.path.join(args.workdir, f)
            circular = True
            break
        if f.startswith('Contigs_1_'):
            draftMito = os.path.join(args.workdir, f)
            break
        if f.startswith('Uncircularized_assemblies_'):
            draftMito = os.path.join(args.workdir, f)
            break
    if circular:
        status('NOVOplasty assembled complete circular genome')
        if args.starting:
            status('Rotating assembly to start with {}'.format(args.starting))
        else:
            status('Rotating assembly to start with Cytochrome b (cob) gene')
        orient_to_start(draftMito, args.out, folder=args.workdir, start=args.starting)
    else:
        numContigs = 0
        contigLength = 0
        with open(args.out, 'w') as outfile:
            with open(draftMito, 'r') as infile:
                for title, seq in SimpleFastaParser(infile):
                    numContigs += 1
                    contigLength += len(seq)
                    outfile.write('>contig_{}\n{}\n'.format(numContigs, softwrap(seq)))
        status('NOVOplasty assembled {} contigs consiting of {:,} bp, but was unable to circularize genome'.format(numContigs, contigLength))

    status('AAFTF mito complete: {}'.format(args.out))
    if not args.pipe:
        shutil.rmtree(args.workdir)

if __name__ == "__main__":
    main()
