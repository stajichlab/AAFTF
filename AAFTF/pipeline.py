"""Support pipelining of AAFTF to simplify all-in-on runs."""

import sys
from argparse import Namespace

import AAFTF.assemble as assemble
import AAFTF.assess as assess
import AAFTF.filter as aaftf_filter
import AAFTF.mito as mito
import AAFTF.polish as polish
import AAFTF.rmdup as rmdup
import AAFTF.sort as aaftf_sort
import AAFTF.sourpurge as sourpurge
import AAFTF.trim as trim
import AAFTF.vecscreen as vecscreen
from AAFTF.utility import checkfile, getRAM, status


def run(parser, args):
    """Script runs entire AAFTF pipeline."""
    # script to run entire AAFTF pipeline
    args_dict = vars(args)
    basename = args_dict['basename']
    RAM = round(0.75*getRAM())
    if not args.memory:
        args_dict['memory'] = str(RAM)

    # run trimming with bbduk
    if not checkfile(basename+'_1P.fastq.gz'):
        trimOpts = ['memory', 'left', 'right', 'basename', 'cpus',
                    'debug', 'minlen']
        trimDict = {k: v for (k, v) in args_dict.items() if k in trimOpts}
        trimDict['method'] = 'bbduk'
        trimDict['pipe'] = True
        trimDict['avgqual'] = 10
        trimargs = Namespace(**trimDict)
        trim.run(parser, trimargs)
    else:
        if args.right:
            status('AAFTF trim output found: {:} {:}'.format(
                basename + '_1P.fastq.gz',
                basename + '_2P.fastq.gz'))
        else:
            status('AAFTF trim output found: {:}'.format(
                basename + '_1P.fastq.gz'))
    if not checkfile(basename + '_1P.fastq.gz'):
        status('AATFT trim failed')
        sys.exit(1)

    # run mitochondrial assembly on bbduk trimmed reads
    if args.right:
        if not checkfile(basename+'.mito.fasta'):
            mitoOpts = ['left', 'right', 'out', 'minlen', 'maxlen', 'seed',
                        'starting', 'workdir', 'pipe', 'reference']
            mitoDict = {k: v for (k, v) in args_dict.items() if k in mitoOpts}
            mitoDict['left'] = basename + '_1P.fastq.gz'
            mitoDict['right'] = basename + '_2P.fastq.gz'
            mitoDict['out'] = basename + '.mito.fasta'
            mitoDict['minlen'] = 10000
            mitoDict['maxlen'] = 100000
            for x in mitoOpts:
                if x not in mitoDict:
                    mitoDict[x] = False
            mitoargs = Namespace(**mitoDict)
            mito.run(parser, mitoargs)
        else:
            status('AAFTF mito output: {}'.format(
                basename + '.mito.fasta'))
    else:
        status('AAFTF mito requires PE reads, ' +
               'skipping mitochondrial de novo assembly')

    # run filtering with bbduk
    if not checkfile(basename+'_filtered_1.fastq.gz'):
        filterOpts = ['screen_accessions', 'screen_urls', 'basename',
                      'cpus', 'debug', 'memory', 'AAFTF_DB', 'workdir']
        filterDict = {k: v for (k, v) in args_dict.items() if k in filterOpts}
        filterDict['aligner'] = 'bbduk'
        filterDict['left'] = basename + '_1P.fastq.gz'
        if args.right:
            filterDict['right'] = basename + '_2P.fastq.gz'
        filterDict['pipe'] = True
        if checkfile(basename + '.mito.fasta'):
            filterDict['screen_local'] = [basename + '.mito.fasta']
        filterargs = Namespace(**filterDict)
        aaftf_filter.run(parser, filterargs)
    else:
        if args.right:
            status('AAFTF filter output found: {:} {:}'.format(
                basename+'_filtered_1.fastq.gz',
                basename+'_filtered_2.fastq.gz'))
        else:
            status('AAFTF filter output found: {:}'.format(
                basename+'_filtered_1.fastq.gz'))

    if not checkfile(basename+'_filtered_1.fastq.gz'):
        status('AATFT filter failed')
        sys.exit(1)

    # run assembly with spades
    if not checkfile(basename+'.spades.fasta'):
        assembleOpts = ['memory', 'cpus', 'debug', 'workdir', 'method',
                        'assembler_args', 'tmpdir']
        asmDict = {k: v for (k, v) in args_dict.items() if k in assembleOpts}
        asmDict['left'] = basename + '_filtered_1.fastq.gz'
        if args.right:
            asmDict['right'] = basename + '_filtered_2.fastq.gz'
        asmDict['out'] = basename + '.spades.fasta'
        asmDict['spades_tmpdir'] = None
        asmDict['pipe'] = True
        asmDict['isolate'] = False
        asmDict['careful'] = True
        asmDict['merged'] = False
        assembleargs = Namespace(**asmDict)
        assemble.run(parser, assembleargs)
    else:
        status('AAFTF assemble output found: {:}'.format(
            basename + '.spades.fasta'))
    if not checkfile(basename + '.spades.fasta'):
        status('AATFT assemble failed')
        sys.exit(1)

    # run vecscreen
    if not checkfile(basename + '.vecscreen.fasta'):
        vecOpts = ['cpus', 'debug', 'workdir', 'AAFTF_DB']
        vecDict = {k: v for (k, v) in args_dict.items() if k in vecOpts}
        vecDict['percent_id'] = False
        vecDict['stringency'] = 'high'
        vecDict['infile'] = basename + '.spades.fasta'
        vecDict['outfile'] = basename + '.vecscreen.fasta'
        vecDict['pipe'] = True
        vecargs = Namespace(**vecDict)
        vecscreen.run(parser, vecargs)
    else:
        status('AAFTF vecscreen output found: {:}'.format(
            basename + '.vecscreen.fasta'))
    if not checkfile(basename + '.vecscreen.fasta'):
        status('AATFT vecscreen failed')
        sys.exit(1)

    # run sourmash purge
    if not checkfile(basename+'.sourpurge.fasta'):
        sourOpts = ['cpus', 'debug', 'workdir', 'AAFTF_DB',
                    'phylum', 'sourdb', 'mincovpct']
        sourDict = {k: v for (k, v) in args_dict.items() if k in sourOpts}
        sourDict['left'] = basename + '_filtered_1.fastq.gz'
        if args.right:
            sourDict['right'] = basename + '_filtered_2.fastq.gz'
        sourDict['input'] = basename + '.vecscreen.fasta'
        sourDict['outfile'] = basename + '.sourpurge.fasta'
        sourDict['kmer'] = 31
        sourDict['taxonomy'] = False
        sourDict['pipe'] = True
        sourDict['sourdb_type'] = 'gbk'
        sourargs = Namespace(**sourDict)
        sourpurge.run(parser, sourargs)
    else:
        status('AAFTF sourpurge output found: {:}'.format(
            basename + '.sourpurge.fasta'))
    if not checkfile(basename + '.sourpurge.fasta'):
        status('AATFT sourpurge failed')
        sys.exit(1)

    # run remove duplicates
    if not checkfile(basename+'.rmdup.fasta'):
        rmdupOpts = ['cpus', 'debug', 'workdir']
        rmdupDict = {k: v for (k, v) in args_dict.items() if k in rmdupOpts}
        rmdupDict['input'] = basename + '.sourpurge.fasta'
        rmdupDict['out'] = basename + '.rmdup.fasta'
        rmdupDict['minlen'] = args_dict['mincontiglen']
        rmdupDict['percent_id'] = 95
        rmdupDict['percent_cov'] = 95
        rmdupDict['exhaustive'] = False
        rmdupDict['pipe'] = True
        rmdupargs = Namespace(**rmdupDict)
        rmdup.run(parser, rmdupargs)
    else:
        status('AAFTF rmdup output found: {:}'.format(
            basename+'.rmdup.fasta'))
    if not checkfile(basename + '.rmdup.fasta'):
        status('AATFT rmdup failed')
        sys.exit(1)

    # run polish to error-correct
    if not checkfile(basename+'.polish.fasta'):
        polishOpts = ['cpus', 'debug', 'workdir', 'iterations', 'memory']
        polishDict = {k: v for (k, v) in args_dict.items() if k in polishOpts}
        polishDict['infile'] = basename + '.rmdup.fasta'
        polishDict['outfile'] = basename + '.polish.fasta'
        polishDict['left'] = basename + '_filtered_1.fastq.gz'
        if args.right:
            polishDict['right'] = basename + '_filtered_2.fastq.gz'
        polishDict['pipe'] = True
        polishargs = Namespace(**polishDict)
        polish.run(parser, polishargs)
    else:
        status('AAFTF polish output found: {:}'.format(
            basename + '.polish.fasta'))
    if not checkfile(basename + '.polish.fasta'):
        status('AATFT polish failed')
        sys.exit(1)

    # sort and rename
    if not checkfile(basename + '.final.fasta'):
        sortDict = {'input': basename + '.polish.fasta',
                    'out':   basename + '.final.fasta',
                    'name':  'scaffold',
                    'minlen': args_dict['mincontiglen']}
        sortargs = Namespace(**sortDict)
        aaftf_sort.run(parser, sortargs)
    else:
        status('AAFTF sort output found: {:}'.format(
            basename + '.final.fasta'))
    if not checkfile(basename + '.final.fasta'):
        status('AATFT sort failed')
        sys.exit(1)

    # assess the assembly
    assessDict = {'input': basename + '.final.fasta',
                  'report': False}
    assessargs = Namespace(**assessDict)
    assess.run(parser, assessargs)
