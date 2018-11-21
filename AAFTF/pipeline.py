import sys
import os
import shutil
import subprocess
from argparse import Namespace 
from AAFTF.utility import status
from AAFTF.utility import printCMD
from AAFTF.utility import getRAM
import AAFTF.trim as trim
import AAFTF.filter as filter
import AAFTF.assemble as assemble
import AAFTF.vecscreen as vecscreen
import AAFTF.sourpurge as sourpurge
import AAFTF.rmdup as rmdup
import AAFTF.pilon as pilon
import AAFTF.sort as sort
import AAFTF.assess as assess

def run(parser,args):
    #script to run entire AAFTF pipeline
    args_dict = vars(args)
    basename = args_dict['basename']
    RAM = round(0.75*getRAM())
    if not args.memory:
        args_dict['memory'] = str(RAM)
    
    #run trimming with bbduk
    if not os.path.isfile(basename+'_1P.fastq.gz'):
        trimOpts = ['memory', 'left', 'right', 'basename', 'cpus', 'debug', 'minlength']
        trimDict = {k:v for (k,v) in args_dict.items() if k in trimOpts}
        trimDict['method'] = 'bbduk'
        trimDict['pipe'] = True
        trimargs = Namespace(**trimDict)
        trim.run(parser, trimargs)
    
    #run filtering with bbduk
    if not os.path.isfile(basename+'_filtered_1.fastq.gz'):
        filterOpts = ['screen_accessions', 'screen_urls', 'basename', 'cpus', 'debug', 'memory', 'AAFTF_DB', 'workdir']
        filterDict = {k:v for (k,v) in args_dict.items() if k in filterOpts}
        filterDict['aligner'] = 'bbduk'
        filterDict['left'] = basename+'_1P.fastq.gz'
        if args.right:
            filterDict['right'] = basename+'_2P.fastq.gz'
        filterDict['pipe'] = True
        filterargs = Namespace(**filterDict)
        filter.run(parser, filterargs)
    
    #run assembly with spades
    if not os.path.isfile(basename+'.spades.fasta'):
        assembleOpts = ['memory', 'cpus', 'debug', 'workdir']
        assembleDict = {k:v for (k,v) in args_dict.items() if k in assembleOpts}
        assembleDict['left'] = basename+'_filtered_1.fastq.gz'
        if args.right:
            assembleDict['right'] = basename+'_filtered_2.fastq.gz'
        assembleDict['out'] = basename+'.spades.fasta'
        assembleDict['spades_tmpdir'] = None
        assembleDict['pipe'] = True
        assembleargs = Namespace(**assembleDict)
        assemble.run(parser, assembleargs)
        
    #run vecscreen
    if not os.path.isfile(basename+'.vecscreen.fasta'):
    	vecOpts = ['cpus', 'debug', 'workdir', 'AAFTF_DB']
    	vecDict = {k:v for (k,v) in args_dict.items() if k in vecOpts}
    	vecDict['percent_id'] = False
    	vecDict['stringency'] = 'high'
    	vecDict['infile'] = basename+'.spades.fasta'
    	vecDict['outfile'] = basename+'.vecscreen.fasta'
    	vecDict['pipe'] = True
    	vecargs = Namespace(**vecDict)
    	vecscreen.run(parser, vecargs)
    
    #run sourmash purge
    if not os.path.isfile(basename+'.sourpurge.fasta'):
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
    
    #run remove duplicates
    if not os.path.isfile(basename+'.rmdup.fasta'):
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
    
    #run pilon to error-correct
    if not os.path.isfile(basename+'.pilon.fasta'):
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
    	
    #sort and rename
    if not os.path.isfile(basename+'.final.fasta'):
    	sortDict = {'input': basename+'.pilon.fasta', 'out': basename+'.final.fasta', 'name': 'scaffold'}
    	sortargs = Namespace(**sortDict)
    	sort.run(parser, sortargs)
    
    #assess the assembly
    assessDict = {'input': basename+'.final.fasta'}
    assessargs = Namespace(**assessDict)
    assess.run(parser, assessargs)
    	
    	
    