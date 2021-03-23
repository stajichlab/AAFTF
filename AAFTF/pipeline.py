import sys
import os
import shutil
import subprocess

from argparse import Namespace 

from AAFTF.utility import status, printCMD, getRAM, checkfile
from AAFTF import trim, assemble, vecscreen, sourcepurge, rmdup, pilon, assess

import AAFTF.filter as aaftf_filter
import AAFTF.sort as aaftf_sort


def run(parser, args):
    #script to run entire AAFTF pipeline
    args_dict = vars(args)
    basename = args_dict['basename']
    RAM = round(0.75*getRAM())
    if not args.memory:
        args_dict['memory'] = str(RAM)
    
    #run trimming with bbduk
    if not checkfile(basename+'_1P.fastq.gz'):
        trimOpts = ['memory', 'left', 'right', 'basename', 'cpus', 'debug', 'minlen']
        trimDict = {k:v for (k,v) in args_dict.items() if k in trimOpts}
        trimDict['method'] = 'bbduk'
        trimDict['pipe'] = True
        trimargs = Namespace(**trimDict)
        trim.run(parser, trimargs)
    else:
    	if args.right:
    	    status('AAFTF trim output found: {:} {:}'.format(basename+'_1P.fastq.gz', basename+'_2P.fastq.gz'))
    	else:
    	    status('AAFTF trim output found: {:}'.format(basename+'_1P.fastq.gz'))
    if not checkfile(basename+'_1P.fastq.gz'):
        status('AATFT trim failed')
        sys.exit(1)
        
    #run filtering with bbduk
    if not checkfile(basename+'_filtered_1.fastq.gz'):
        filterOpts = ['screen_accessions', 'screen_urls', 'basename', 'cpus', 'debug', 'memory', 'AAFTF_DB', 'workdir']
        filterDict = {k:v for (k,v) in args_dict.items() if k in filterOpts}
        filterDict['aligner'] = 'bbduk'
        filterDict['left'] = basename+'_1P.fastq.gz'
        if args.right:
            filterDict['right'] = basename+'_2P.fastq.gz'
        filterDict['pipe'] = True
        filterargs = Namespace(**filterDict)
        aaftf_filter.run(parser, filterargs)
    else:
    	if args.right:
            status('AAFTF filter output found: {:} {:}'.format(basename+'_filtered_1.fastq.gz', basename+'_filtered_2.fastq.gz'))
    	else:
            status('AAFTF filter output found: {:}'.format(basename+'_filtered_1.fastq.gz'))
    if not checkfile(basename+'_filtered_1.fastq.gz'):
        status('AATFT filter failed')
        sys.exit(1)
        
    #run assembly with spades
    if not checkfile(basename+'.spades.fasta'):
        assembleOpts = ['memory', 'cpus', 'debug', 'workdir', 'method', 'assembler_args', 'tmpdir']
        assembleDict = {k:v for (k,v) in args_dict.items() if k in assembleOpts}
        assembleDict['left'] = basename+'_filtered_1.fastq.gz'
        if args.right:
            assembleDict['right'] = basename+'_filtered_2.fastq.gz'
        assembleDict['out'] = basename+'.spades.fasta'
        assembleDict['spades_tmpdir'] = None
        assembleDict['pipe'] = True
        assembleargs = Namespace(**assembleDict)
        assemble.run(parser, assembleargs)
    else:
    	status('AAFTF assemble output found: {:}'.format(basename+'.spades.fasta'))
    if not checkfile(basename+'.spades.fasta'):
        status('AATFT assemble failed')
        sys.exit(1)
        
    #run vecscreen
    if not checkfile(basename+'.vecscreen.fasta'):
        vecOpts = ['cpus', 'debug', 'workdir', 'AAFTF_DB']
        vecDict = {k:v for (k,v) in args_dict.items() if k in vecOpts}
        vecDict['percent_id'] = False
        vecDict['stringency'] = 'high'
        vecDict['infile'] = basename+'.spades.fasta'
        vecDict['outfile'] = basename+'.vecscreen.fasta'
        vecDict['pipe'] = True
        vecargs = Namespace(**vecDict)
        vecscreen.run(parser, vecargs)
    else:
    	status('AAFTF vecscreen output found: {:}'.format(basename+'.vecscreen.fasta'))
    if not checkfile(basename+'.vecscreen.fasta'):
        status('AATFT vecscreen failed')
        sys.exit(1)
    
    #run sourmash purge
    if not checkfile(basename+'.sourpurge.fasta'):
        sourOpts = ['cpus', 'debug', 'workdir', 'AAFTF_DB', 'phylum', 'sourdb', 'mincovpct']
        sourDict = {k:v for (k,v) in args_dict.items() if k in sourOpts}
        sourDict['left'] = basename+'_filtered_1.fastq.gz'
        if args.right:
            sourDict['right'] = basename+'_filtered_2.fastq.gz'
        sourDict['input'] = basename+'.vecscreen.fasta'
        sourDict['outfile'] = basename+'.sourpurge.fasta'
        sourDict['taxonomy'] = False
        sourDict['pipe'] = True
        sourargs = Namespace(**sourDict)
        sourpurge.run(parser, sourargs)
    else:
    	status('AAFTF sourpurge output found: {:}'.format(basename+'.sourpurge.fasta'))
    if not checkfile(basename+'.sourpurge.fasta'):
        status('AATFT sourpurge failed')
        sys.exit(1)

    #run remove duplicates
    if not checkfile(basename+'.rmdup.fasta'):
        rmdupOpts = ['cpus', 'debug', 'workdir']
        rmdupDict = {k:v for (k,v) in args_dict.items() if k in rmdupOpts}
        rmdupDict['input'] = basename+'.sourpurge.fasta'
        rmdupDict['out'] = basename+'.rmdup.fasta'
        rmdupDict['minlen'] = args_dict['mincontiglen']
        rmdupDict['percent_id'] = 95
        rmdupDict['percent_cov'] = 95
        rmdupDict['exhaustive'] = False
        rmdupDict['pipe'] = True
        rmdupargs = Namespace(**rmdupDict)
        rmdup.run(parser, rmdupargs)
    else:
    	status('AAFTF rmdup output found: {:}'.format(basename+'.rmdup.fasta'))
    if not checkfile(basename+'.rmdup.fasta'):
        status('AATFT rmdup failed')
        sys.exit(1)
            
    #run pilon to error-correct
    if not checkfile(basename+'.pilon.fasta'):
        pilonOpts = ['cpus', 'debug', 'workdir', 'iterations']
        pilonDict = {k:v for (k,v) in args_dict.items() if k in pilonOpts}
        pilonDict['infile'] = basename+'.rmdup.fasta'
        pilonDict['outfile'] = basename+'.pilon.fasta'
        pilonDict['left'] = basename+'_filtered_1.fastq.gz'
        if args.right:
            pilonDict['right'] = basename+'_filtered_2.fastq.gz'
        pilonDict['pipe'] = True
        pilonargs = Namespace(**pilonDict)
        pilon.run(parser, pilonargs)
    else:
    	status('AAFTF pilon output found: {:}'.format(basename+'.pilon.fasta'))
    if not checkfile(basename+'.pilon.fasta'):
        status('AATFT pilon failed')
        sys.exit(1)
        
    #sort and rename
    if not checkfile(basename+'.final.fasta'):
        sortDict = {'input': basename+'.pilon.fasta', 'out': basename+'.final.fasta', 'name': 'scaffold'}
        sortargs = Namespace(**sortDict)
        aaftf_sort.run(parser, sortargs)
    else:
    	status('AAFTF sort output found: {:}'.format(basename+'.final.fasta'))
    if not checkfile(basename+'.final.fasta'):
        status('AATFT sort failed')
        sys.exit(1)
   
    #assess the assembly
    assessDict = {'input': basename+'.final.fasta', 'report': False}
    assessargs = Namespace(**assessDict)
    assess.run(parser, assessargs)
        
        
    
