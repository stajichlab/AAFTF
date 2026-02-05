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

    # Helper function to create namespace with required defaults
    def create_namespace(options, required_args=None, **extra_args):
        """Create a Namespace with filtered options and required arguments."""
        namespace_dict = {k: v for (k, v) in args_dict.items() if k in options}
        if required_args:
            namespace_dict.update(required_args)
        namespace_dict.update(extra_args)
        return Namespace(**namespace_dict)

    # Helper function to check step output and handle failures
    def check_step_success(output_file, step_name):
        """Check if step completed successfully."""
        if not checkfile(output_file):
            status(f'AAFTF {step_name} failed')
            sys.exit(1)
        return True

    # run trimming with bbduk
    if not checkfile(basename+'_1P.fastq.gz'):
        trimOpts = ['memory', 'left', 'right', 'basename', 'cpus',
                    'debug', 'minlen']
        trim_args = create_namespace(
            trimOpts,
            required_args={'method': 'bbduk', 'pipe': True, 'avgqual': 10}
        )
        trim.run(parser, trim_args)
    else:
        if args.right:
            status('AAFTF trim output found: {:} {:}'.format(
                basename + '_1P.fastq.gz',
                basename + '_2P.fastq.gz'))
        else:
            status('AAFTF trim output found: {:}'.format(
                basename + '_1P.fastq.gz'))
    check_step_success(basename + '_1P.fastq.gz', 'trim')

    # run mitochondrial assembly on bbduk trimmed reads
    if args.right:
        if not checkfile(basename+'.mito.fasta'):
            mitoOpts = ['left', 'right', 'out', 'minlen', 'maxlen', 'seed',
                        'starting', 'workdir', 'pipe', 'reference']
            mito_args = create_namespace(
                mitoOpts,
                required_args={
                    'left': basename + '_1P.fastq.gz',
                    'right': basename + '_2P.fastq.gz',
                    'out': basename + '.mito.fasta',
                    'minlen': 10000,
                    'maxlen': 100000,
                    'pipe': True
                },
                **{x: False for x in mitoOpts if x not in args_dict}
            )
            mito.run(parser, mito_args)
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
        filter_args = create_namespace(
            filterOpts,
            required_args={
                'aligner': 'bbduk',
                'left': basename + '_1P.fastq.gz',
                'pipe': True
            }
        )
        if args.right:
            filter_args.right = basename + '_2P.fastq.gz'
        if checkfile(basename + '.mito.fasta'):
            filter_args.screen_local = [basename + '.mito.fasta']
        aaftf_filter.run(parser, filter_args)
    else:
        if args.right:
            status('AAFTF filter output found: {:} {:}'.format(
                basename+'_filtered_1.fastq.gz',
                basename+'_filtered_2.fastq.gz'))
        else:
            status('AAFTF filter output found: {:}'.format(
                basename+'_filtered_1.fastq.gz'))
    check_step_success(basename+'_filtered_1.fastq.gz', 'filter')

    # run assembly with specified method
    assembly_method = args.method if hasattr(args, 'method') else 'spades'
    assembly_file = basename + f'.{assembly_method}.fasta'
    if not checkfile(assembly_file):
        assembleOpts = ['memory', 'cpus', 'debug', 'workdir', 'method',
                        'assembler_args', 'tmpdir']
        asm_args = create_namespace(
            assembleOpts,
            required_args={
                'left': basename + '_filtered_1.fastq.gz',
                'out': assembly_file,
                'pipe': True,
                'method': assembly_method
            }
        )
        if args.right:
            asm_args.right = basename + '_filtered_2.fastq.gz'
        # Set assembly-specific parameters
        if assembly_method == 'spades':
            asm_args.spades_tmpdir = None
            asm_args.isolate = False
            asm_args.careful = True
        asm_args.merged = False
        assemble.run(parser, asm_args)
    else:
        status('AAFTF assemble output found: {:}'.format(assembly_file))
    check_step_success(assembly_file, 'assemble')

    # run vecscreen
    vecscreen_file = basename + '.vecscreen.fasta'
    if not checkfile(vecscreen_file):
        vecOpts = ['cpus', 'debug', 'workdir', 'AAFTF_DB']
        vec_args = create_namespace(
            vecOpts,
            required_args={
                'percent_id': False,
                'stringency': 'high',
                'infile': assembly_file,
                'outfile': vecscreen_file,
                'pipe': True
            }
        )
        vecscreen.run(parser, vec_args)
    else:
        status('AAFTF vecscreen output found: {:}'.format(vecscreen_file))
    check_step_success(vecscreen_file, 'vecscreen')

    # run sourmash purge
    sourpurge_file = basename + '.sourpurge.fasta'
    if not checkfile(sourpurge_file):
        sourOpts = ['cpus', 'debug', 'workdir', 'AAFTF_DB',
                    'phylum', 'sourdb', 'mincovpct']
        sour_args = create_namespace(
            sourOpts,
            required_args={
                'left': basename + '_filtered_1.fastq.gz',
                'input': vecscreen_file,
                'outfile': sourpurge_file,
                'kmer': 31,
                'taxonomy': False,
                'pipe': True,
                'sourdb_type': 'gbk'
            }
        )
        if args.right:
            sour_args.right = basename + '_filtered_2.fastq.gz'
        sourpurge.run(parser, sour_args)
    else:
        status('AAFTF sourpurge output found: {:}'.format(sourpurge_file))
    check_step_success(sourpurge_file, 'sourpurge')

    # run remove duplicates
    rmdup_file = basename+'.rmdup.fasta'
    if not checkfile(rmdup_file):
        rmdupOpts = ['cpus', 'debug', 'workdir']
        rmdup_args = create_namespace(
            rmdupOpts,
            required_args={
                'input': sourpurge_file,
                'out': rmdup_file,
                'minlen': args.mincontiglen,
                'percent_id': 95,
                'percent_cov': 95,
                'exhaustive': False,
                'pipe': True
            }
        )
        rmdup.run(parser, rmdup_args)
    else:
        status('AAFTF rmdup output found: {:}'.format(rmdup_file))
    check_step_success(rmdup_file, 'rmdup')

    # run polish to error-correct
    polish_file = basename+'.polish.fasta'
    if not checkfile(polish_file):
        polishOpts = ['cpus', 'debug', 'workdir', 'iterations', 'memory']
        polish_args = create_namespace(
            polishOpts,
            required_args={
                'infile': rmdup_file,
                'outfile': polish_file,
                'left': basename + '_filtered_1.fastq.gz',
                'pipe': True
            }
        )
        if args.right:
            polish_args.right = basename + '_filtered_2.fastq.gz'
        polish.run(parser, polish_args)
    else:
        status('AAFTF polish output found: {:}'.format(polish_file))
    check_step_success(polish_file, 'polish')

    # sort and rename
    final_file = basename + '.final.fasta'
    if not checkfile(final_file):
        sort_args = Namespace(
            input=polish_file,
            out=final_file,
            name='scaffold',
            minlen=args.mincontiglen
        )
        aaftf_sort.run(parser, sort_args)
    else:
        status('AAFTF sort output found: {:}'.format(final_file))
    check_step_success(final_file, 'sort')

    # assess the assembly
    assess_args = Namespace(
        input=final_file,
        report=False
    )
    assess.run(parser, assess_args)