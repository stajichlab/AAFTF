"""Run the Mitochondria assembly tool NOVOPlasty."""

import os
import shutil
import subprocess
import sys
import uuid

from Bio.SeqIO.FastaIO import SimpleFastaParser

from AAFTF.resources import Mitoseqs
from AAFTF.utility import (GuessRL, RevComp, execute, getRAM, printCMD,
                           softwrap, status, which_path)


def orient_to_start(fasta_in, fasta_out, folder='.', start=False):
    """Reorient the MT assembly based on a starting gene (if found)."""
    # if not starting, then use cytochrome oxidase (cob)
    startFile = os.path.join(folder, f'{uuid.uuid4()}.fasta')
    if not start:
        # generated as spoa consensus from select fungal cob genes
        # move this to a configurable file
        cob1 = Mitoseqs['COB1']
        with open(startFile, 'w') as outfile:
            outfile.write(f'>COB\n{softwrap(cob1)}\n')
    else:
        shutil.copyfile(start, startFile)

    # load sequence into dictionary
    initial_seq = ''
    # header = ''
    with open(fasta_in) as infile:
        for title, seq in SimpleFastaParser(infile):
            initial_seq = seq
            # header = title

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
        status('ERROR: unable to rotate because did ' +
               'not find --starting sequence\n')
        with open(fasta_out, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format('mt', softwrap(initial_seq)))
    elif len(alignments) > 1:
        status('ERROR: unable to rotate because found multiple alignments\n')
        for x in alignments:
            sys.stderr.write(f'{x}\n')
        with open(fasta_out, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format('mt', softwrap(initial_seq)))
    if os.path.isfile(startFile):
        os.remove(startFile)


def run(parser, args):
    """Run the NOVOplasty tool."""
    # first check if NOVOplasty and minimap2 are installed, else exit
    programs = ['NOVOPlasty.pl', 'minimap2']
    for x in programs:
        if not which_path(x):
            status(f'ERROR: {x} is not installed, exiting')
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
            seedFasta = os.path.abspath(os.path.join(
                os.path.dirname(__file__), 'data', 'mito-seed.fasta'))
        else:
            seedFasta = os.path.abspath(args.reference)
    else:
        seedFasta = os.path.abspath(args.seed)

    # now write the novoplasty config file
    defaultConfig = os.path.join(os.path.dirname(
        __file__), 'data', 'novoplasty-config.txt')
    novoConfig = os.path.join(args.workdir, 'novo-config.txt')
    if args.reference:
        refgenome = os.path.abspath(args.reference)
    else:
        refgenome = ''
    checkWords = ("<PROJECT>", "<MINLEN>", "<MAXLEN>", "<MAXMEM>",
                  "<SEED>", "<READLEN>", "<FORWARD>", "<REVERSE>",
                  "<REFERENCE>")
    repWords = (unique_id,               # project
                str(args.minlen),        # minlen
                str(args.maxlen),        # maxlen
                str(int(getRAM()*.75)),  # maxRAM
                seedFasta,               # seed fasta seq
                str(read_len),           # read length
                os.path.abspath(args.left),   # forward read
                os.path.abspath(args.right),  # rev read
                refgenome)               # ref genome file
    with open(novoConfig, 'w') as outfile:
        with open(defaultConfig) as infile:
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
        p1 = subprocess.Popen(cmd, cwd=args.workdir,
                              stdout=logfile, stderr=logfile)
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
            status(f'Rotating assembly to start with {args.starting}')
        else:
            status('Rotating assembly to start with Cytochrome b (cob) gene')
        orient_to_start(draftMito, args.out,
                        folder=args.workdir, start=args.starting)
    else:
        numContigs = 0
        contigLength = 0
        with open(args.out, 'w') as outfile:
            with open(draftMito) as infile:
                for title, seq in SimpleFastaParser(infile):
                    numContigs += 1
                    contigLength += len(seq)
                    outfile.write('>contig_{}\n{}\n'.format(
                        numContigs, softwrap(seq)))
        # weird formatting here for PEP8
        status(
            'NOVOplasty assembled {} contigs consisting of {:,} bp,'.format(
                numContigs, contigLength) +
            'but was unable to circularize genome')

    status(f'AAFTF mito complete: {args.out}')
    if not args.pipe:
        shutil.rmtree(args.workdir)
