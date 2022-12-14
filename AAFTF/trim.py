"""Trims FASTQ files for reads.

This is usually for Illumina reads to quality trim reads.
This uses either fastp, includes merging step for paired reads, OR
trimmomatic. Expects adaptor sequence files to be in trimmomatic installed folder.
"""

import os
import subprocess
import sys
from os.path import dirname

from AAFTF.utility import (Fzip_inplace, SafeRemove, countfastq, getRAM,
                           printCMD, status, which_path)

TRIMMOMATIC_TRUSEQSE = "adapters/TruSeq3-SE.fa"
TRIMMOMATIC_TRUSEQPE = "adapters/TruSeq3-PE.fa"
# process trimming reads with trimmomatic
# Homebrew install of trimmomatic uses a shell script
'''
#!/bin/bash -l
TRIMJAR=/usr/local/Cellar/trimmomatic/0.36/libexec/trimmomatic-0.36.jar
exec java -jar $TRIMJAR "$@"
'''
# while bioconda install uses a python script that launches java apps


def find_trimmomatic():
    """Finds the trimmomatic jar file."""
    trim_path = which_path('trimmomatic')
    if trim_path:
        with open(os.path.abspath(trim_path)) as trim_shell:
            firstLine = trim_shell.readline()
            if '#!/bin/bash' in firstLine:  # homebrew get jar location
                for line in trim_shell:
                    if line.startswith('exec java'):
                        items = line.split(' ')
                        for x in items:
                            if x.endswith('.jar'):
                                return x
            elif '#!/usr/bin/env python' in firstLine:
                trimjardir = os.path.dirname(os.path.realpath(trim_path))
                return os.path.join(trimjardir, 'trimmomatic.jar')
            else:
                return False
    else:
        return False

# flake8: noqa: C901
def run(parser, args):
    """Run command for the module subtool of AAFTF."""
    if not args.basename:
        if '_' in os.path.basename(args.left):
            args.basename = os.path.basename(args.left).split('_')[0]
        elif '.' in os.path.basename(args.left):
            args.basename = os.path.basename(args.left).split('.')[0]
        else:
            args.basename = os.path.basename(args.left)

    total = countfastq(args.left)
    if args.right:
        total = total*2
    status(f'Loading {total:,} total reads')

    DEVNULL = open(os.devnull, 'w')
    if args.method == 'bbduk':
        if args.memory:
            MEM = f'-Xmx{args.memory}g'
        else:
            MEM = f'-Xmx{round(0.6*getRAM())}g'

        status('Adapter trimming using BBDuk')
        cmd = ['bbduk.sh', MEM,
               'ref=adapters',
               f't={args.cpus}',
               'ktrim=r',
               'k=23',
               'mink=11',
               f'minlen={args.minlen}',
               'hdist=1',
               f'maq={args.avgqual}',
               'ftm=5',
               'tpe',
               'tbo',
               'overwrite=true']
        if args.left and args.right:
            cmd += [f'in1={args.left}',
                    f'in2={args.right}',
                    f'out1={args.basename}_1P.fastq.gz',
                    f'out2={args.basename}_2P.fastq.gz']
        elif args.left:
            cmd += [f'in={args.left}',
                    f'out={args.basename}_1U.fastq.gz']

        printCMD(cmd)
        if args.debug:
            subprocess.run(cmd)
        else:
            subprocess.run(cmd, stderr=DEVNULL)

        if args.right:
            clean = countfastq(f'{args.basename}_1P.fastq.gz')
            clean = clean*2
            status(f'{clean:,} reads remaining and writing to file')
            status('Trimming finished:\n\tFor: {:}\n\tRev {:}'.format(
                args.basename + '_1P.fastq.gz',
                args.basename + '_2P.fastq.gz'))
            if not args.pipe:
                status('Your next command might be:\n\t' +
                       'AAFTF filter -l {:} -r {:} -o {:} -c {:}\n'.format(
                           args.basename+'_1P.fastq.gz',
                           args.basename+'_2P.fastq.gz',
                           args.basename,
                           args.cpus))
        else:
            clean = countfastq(f'{args.basename}_1U.fastq.gz')
            status(f'{clean:,} reads remaining and writing to file')
            status('Trimming finished:\n\tSingle: {:}'.format(
                args.basename+'_1U.fastq.gz'))
            if not args.pipe:
                status('Your next command might be:\n\t' +
                       'AAFTF filter -l {:} -o {:} -c {:}\n'.format(
                           args.basename+'_1U.fastq.gz',
                           args.basename,
                           args.cpus))

    elif args.method == 'trimmomatic':
        # find path
        trimmomatic_path = find_trimmomatic()
        if trimmomatic_path:
            jarfile = trimmomatic_path
        elif args.trimmomatic:
            jarfile = args.trimmomatic
        else:
            status('Trimmomatic cannot be found - ' +
                   'please provide location of trimmomatic.jar file.')
            sys.exit(1)

        if jarfile:
            path_to_adaptors = args.trimmomatic_adaptors
            leadingwindow = "LEADING:%d" % (args.trimmomatic_leadingwindow)
            trailingwindow = "TRAILING:%d" % (args.trimmomatic_trailingwindow)
            slidingwindow = "SLIDINGWINDOW:%s" % (
                args.trimmomatic_slidingwindow)

            quality = args.trimmomatic_quality
            quality = "-%s" % (quality)  # add leading dash

            if not os.path.exists(path_to_adaptors):
                if args.right:
                    path_to_adaptors = os.path.join(dirname(jarfile),
                                                    TRIMMOMATIC_TRUSEQPE)
                else:
                    path_to_adaptors = os.path.join(dirname(jarfile),
                                                    TRIMMOMATIC_TRUSEQSE)

                if not os.path.exists(path_to_adaptors):
                    findpath = dirname(jarfile)
                    path_to_adaptors = ""
                    while findpath:
                        if os.path.exists(findpath + "/share"):
                            if args.right:
                                path_to_adaptors = os.path.join(
                                    findpath,
                                    "/share/trimmomatic",
                                    TRIMMOMATIC_TRUSEQPE)
                            else:
                                path_to_adaptors = os.path.join(
                                    findpath,
                                    "/share/trimmomatic",
                                    TRIMMOMATIC_TRUSEQSE)
                            break
                        findpath = dirname(findpath)

                if not os.path.exists(path_to_adaptors):
                    status("Cannot find adaptors file please specify manually")
                    return

            clipstr = args.trimmomatic_clip % (path_to_adaptors)

            cmd = []

            if args.left and args.right:
                cmd = ['java', '-jar', jarfile, 'PE',
                       '-threads', str(args.cpus), quality,
                       args.left, args.right,
                       args.basename+'_1P.fastq',
                       args.basename+'_1U.fastq',
                       args.basename+'_2P.fastq',
                       args.basename+'_2U.fastq',
                       clipstr, leadingwindow, trailingwindow, slidingwindow,
                       "MINLEN:%d" % (args.minlen)]
            elif args.left and not args.right:
                cmd = ['java', '-jar', jarfile, 'SE',
                       '-threads', str(args.cpus),
                       quality,  args.left,
                       args.basename+'_1U.fastq',
                       clipstr, leadingwindow, trailingwindow, slidingwindow,
                       "MINLEN:%d" % (args.minlen)]
            else:
                status("Must provide left and right pairs or single read set")
                return

            status('Running trimmomatic adapter and quality trimming')
            printCMD(cmd)
            if args.debug:
                subprocess.run(cmd)
            else:
                subprocess.run(cmd, stderr=DEVNULL)
            if args.right:
                status('Compressing trimmed PE FASTQ files')
                Fzip_inplace(args.basename+'_1P.fastq', args.cpus)
                Fzip_inplace(args.basename+'_2P.fastq', args.cpus)
                SafeRemove(args.basename+'_1U.fastq')
                SafeRemove(args.basename+'_2U.fastq')
                status('Trimming finished:\n\tFor: {:}\n\tRev {:}'.format(
                    args.basename+'_1P.fastq.gz',
                    args.basename+'_2P.fastq.gz'))
                if not args.pipe:
                    status('Your next command might be:\n\t' +
                           'AAFTF filter -l {:} -r {:} -o {:} -c {:}\n'.format(
                               args.basename+'_1P.fastq.gz',
                               args.basename+'_2P.fastq.gz',
                               args.basename,
                               args.cpus))
            else:
                status('Compressing trimmed SE FASTQ file')
                Fzip_inplace(args.basename + '_1U.fastq', args.cpus)
                status('Trimming finished:\n\tSingle: {:}'.format(
                    args.basename + '_1U.fastq.gz'))
                if not args.pipe:
                    status('Your next command might be:\n\t' +
                           'AAFTF filter -l {:} -o {:} -c {:}\n'.format(
                               args.basename+'_1U.fastq.gz',
                               args.basename,
                               args.cpus))

    elif args.method == 'fastp':
        status('Adapter trimming using fastp')
        cmd = ['fastp', '--low_complexity_filter',
               '-l', f'{args.minlen}',
               '--average_qual', f'{args.avgqual}',
               '-w', f'{args.cpus}']

#               '-wref=adapters', 't={:}'.format(args.cpus), 'ktrim=r',
#           'k=23', 'mink=11', 'minlen={:}'.format(args.minlen), 'hdist=1',
#           'ftm=5', 'tpe', 'tbo', 'overwrite=true']
        if args.left and args.right:
            # could add merging ...
            cmd += [f'--in1={args.left}',
                    f'--in2={args.right}',
                    f'--out1={args.basename}_1P.fastq.gz',
                    f'--out2={args.basename}_2P.fastq.gz'
                    ]
            if args.merge:
                cmd += ['--merge',
                        f'--merged_out={args.basename}_MG.fastq.gz']

        elif args.left:
            cmd += [f'--in={args.left}',
                    f'--out={args.basename}_1U.fastq.gz']
        if args.dedup:
            cmd += ['--dedup']
        if args.cutfront:
            cmd += ['--cut_front']
        if args.cuttail:
            cmd += ['--cut_tail']
        if args.cutright:
            cmd += ['--cut_right']

        cmd += [f'--html={args.basename}.fastp.html',
                f'--json={args.basename}.fastp.json']
        printCMD(cmd)
        if args.debug:
            subprocess.run(cmd)
        else:
            subprocess.run(cmd, stderr=DEVNULL)

        if args.right:
            clean = countfastq(f'{args.basename}_1P.fastq.gz')
            clean = clean*2
            status(f'{clean:,} reads remaining and writing to file')
            status('Trimming finished:\n\tFor: {:}\n\tRev {:}'.format(
                args.basename+'_1P.fastq.gz',
                args.basename+'_2P.fastq.gz'))
            if not args.pipe:
                status('Your next command might be:\n\t' +
                       'AAFTF filter -l {:} -r {:} -o {:} -c {:}\n'.format(
                           args.basename+'_1P.fastq.gz',
                           args.basename+'_2P.fastq.gz',
                           args.basename,
                           args.cpus))
        else:
            clean = countfastq(f'{args.basename}_1U.fastq.gz')
            status(f'{clean:,} reads remaining and writing to file')
            status('Trimming finished:\n\tSingle: {:}'.format(
                args.basename + '_1U.fastq.gz'))
            if not args.pipe:
                status('Your next command might be:\n\t' +
                       'AAFTF filter --left {:} -o {:} -c {:}\n'.format(
                           args.basename+'_1U.fastq.gz',
                           args.basename, args.cpus))

    else:
        status(f'Uknown trimming method: {args.method}')
