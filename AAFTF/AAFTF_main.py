#!/usr/bin/env python3
"""AAFTF main framework for submodule run."""

# note structure of code taken from poretools
# https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import argparse as ap
import sys

# AAFTF imports
from AAFTF.__version__ import __version__  # noqa: E402
from AAFTF.utility import status

myversion = __version__


def run_subtool(parser, args):
    """Run the subtool of the AAFTF pipeline."""
    if args.command == 'runall':
        print("runall")
    elif args.command == 'trim':
        import AAFTF.trim as submodule
    elif args.command == 'filter':
        import AAFTF.filter as submodule
    elif args.command == 'assemble':
        import AAFTF.assemble as submodule
    elif args.command == 'vecscreen':
        import AAFTF.vecscreen as submodule
    elif args.command == 'fcs_screen':
        import AAFTF.fcs_screen as submodule
    elif args.command == 'fcs_gx_purge':
        import AAFTF.fcs_gx_purge as submodule
    elif args.command == 'sourpurge':
        import AAFTF.sourpurge as submodule
    elif args.command == 'rmdup':
        import AAFTF.rmdup as submodule
    elif args.command == 'polish':
        import AAFTF.polish as submodule
    elif args.command == 'pilon':
        import AAFTF.polish as submodule
    elif args.command == 'assess':
        import AAFTF.assess as submodule
    elif args.command == 'sort':
        import AAFTF.sort as submodule
    elif args.command == 'pipeline':
        import AAFTF.pipeline as submodule
    elif args.command == 'mito':
        import AAFTF.mito as submodule
    else:
        parser.parse_args('')
        return
    # run the chosen submodule.
    submodule.run(parser, args)


class ArgumentParserWithDefaults(ap.ArgumentParser):
    """Argument parsing top level."""
    def __init__(self, *args, **kwargs):
        """Init function for ArgumentParser."""
        super().__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet",
                          help="Do not output warnings to stderr",
                          action="store_true",
                          dest="quiet")

    def error(self, message):
        """Error message passing from argument passing."""
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    """Present the main AAFTF module submenus."""
    #########################################
    # create the top-level parser
    #########################################
    parser = ap.ArgumentParser(
        prog='AAFTF',
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q", "--quiet",
                        help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")
    parser.add_argument("-v", "--version", help="Installed AAFTF version",
                        action="version",
                        version="%(prog)s " + str(myversion))

    subparsers = parser.add_subparsers(title='[sub-commands]',
                                       dest='command',
                                       parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################

    ##########
    # trim
    ##########
    # arguments
    # --method bbduk*, trimmomatic, fastp
    # --trimmomatic: arguments are path to JAR or application respectively
    # assume java is PATH already for trimmomatic
    # -o / --outdir: write outdir
    # -p / --prefix: outfile prefix
    # -ml / --minlen: min read length
    # --avgqual (default 10)
    # read info, either paired data are required or singleton
    # --left: left or forward reads
    # --right: right or reverse reads
    #
    # --merge: merge reads in fastp

    parser_trim = subparsers.add_parser(
        'trim',
        aliases=['trim_reads', 'read_trim'],
        description="This command trims reads in FASTQ format to " +
        "remove low quality reads and trim adaptor sequences",
        help='Trim FASTQ input reads')

    parser_trim.add_argument(
        '-o', '--out', type=str,
        required=False, dest='basename',
        help="Output basename, default to base name of --left reads")

    parser_trim.add_argument(
        '-c', '--cpus', type=int,
        metavar="cpus",
        required=False,
        default=1,
        help="Number of CPUs/threads to use.")

    parser_trim.add_argument(
        '-ml', '--minlen', type=int,
        default=75,
        required=False,
        help="Minimum read length after trimming, default: 75")

    parser_trim.add_argument(
        '-l', '--left', type=str,
        required=True,
        help='left/forward reads of paired-end FASTQ or single-end FASTQ.')

    parser_trim.add_argument('-r', '--right', type=str,
                             required=False,
                             help='right/reverse reads of paired-end FASTQ.')

    parser_trim.add_argument(
        '-aq', '--avgqual', type=int,
        default=10,
        required=False,
        help="Average Quality of reads must be > than this, default: 10")

    parser_trim.add_argument(
        '--dedup', action='store_true',
        required=False,
        help="Run fastp deuplication of fastq reads (default uses ~4gb mem)" +
        " default: false")

    parser_trim.add_argument(
        '--cutfront', action='store_true',
        required=False,
        help="Run fastp 5' trimming based on quality. " +
        "default: false. \n" +
        "WARNING: this operation will interfere deduplication for SE data")

    parser_trim.add_argument(
        '--cuttail', action='store_true',
        required=False,
        help="Run fastp 3' trimming based on quality. " +
        "default: false.\n" +
        "WARNING: this operation will interfere deduplication for SE data")

    parser_trim.add_argument(
        '--cutright', action='store_true',
        required=False,
        help="Run fastp move a sliding window from front to tail, " +
        "if meet one window with mean quality < threshold. \n" +
        "default: false.\n" +
        "WARNING: this operation will interfere deduplication for SE data")

    parser_trim.add_argument('-v', '--debug', action='store_true',
                             help="Provide debugging messages")

    parser_trim.add_argument('--pipe', action='store_true',
                             help="AAFTF is running in pipeline mode")

    parser_trim.add_argument('--method', default='bbduk',
                             choices=['bbduk', 'trimmomatic', 'fastp'],
                             help='Program to use for adapter trimming')

    parser_trim.add_argument('-m', '--memory', type=int,
                             dest='memory', required=False,
                             help="Max Memory (in GB)")

    parser_trim.add_argument('--merge', action='store_true',
                             help="Merge paired end reads when running fastp")

    parser_trim.add_argument('--prefix', type=str, required=False,
                             help="Prefix for outfiles")

    tool_group = parser_trim.add_mutually_exclusive_group(required=False)

    tool_group.add_argument('--trimmomatic', '--jar',
                            metavar='trimmomatic_jar',
                            type=str, required=False,
                            help='Trimmomatic JAR path')
    trimmomatic_group = parser_trim.add_argument_group(
        title='Trimmomatic options',
        description="Trimmomatic trimming options")

    trimmomatic_group.add_argument(
        '--trimmomatic_adaptors',
        default="TruSeq3-PE.fa",
        help="Trimmomatic adaptor file, default: TruSeq3-PE.fa")

    trimmomatic_group.add_argument(
        '--trimmomatic_clip',
        default="ILLUMINACLIP:%s:2:30:10",
        help="Trimmomatic clipping, " +
        "default: ILLUMINACLIP:TruSeq3-PE.fa:2:30:10")

    trimmomatic_group.add_argument(
        '--trimmomatic_leadingwindow',
        default="3", type=int,
        help="Trimmomatic window processing arguments, default: LEADING:3")

    trimmomatic_group.add_argument(
        '--trimmomatic_trailingwindow',
        default="3", type=int,
        help="Trimmomatic window processing arguments, default: TRAILING:3")

    trimmomatic_group.add_argument(
        '--trimmomatic_slidingwindow',
        default="4:15", type=str,
        help="Trimmomatic window processing arguments, " +
        "default: SLIDINGWINDOW:4:15")

    trimmomatic_group.add_argument(
        '--trimmomatic_quality',
        default="phred33",
        help="Trimmomatic quality encoding -phred33 or phred64")

    ##########
    # mito-asm assembly mitochondrial genome
    ##########

    parser_mito = subparsers.add_parser(
        'mito',
        aliases=['mito_asm', 'mitochondria'],
        description="De novo assembly of mitochondrial genome using " +
        "NOVOplasty, takes PE Illumina adapter trimmed data.",
        help='De novo assembly of mitochondrial genome')
    parser_mito.add_argument('-l', '--left', required=True,
                             help="Left (Forward) reads")

    parser_mito.add_argument('-r', '--right', required=True,
                             help="Right (Reverse) reads")

    parser_mito.add_argument('-o', '--out', type=str,
                             required=True,
                             help="Output FASTA file for mitochondrial genome")

    parser_mito.add_argument('--minlen', default=10000, type=int,
                             help="Minimum expected genome size")

    parser_mito.add_argument('--maxlen', default=100000, type=int,
                             help="Maximum expected genome size")

    parser_mito.add_argument(
        '-s', '--seed', required=False,
        help="Seed sequence, ie related mitochondrial genome." +
        "default: A. nidulans")

    parser_mito.add_argument(
        '--starting', required=False,
        help="FASTA file of start sequence, rotate genome to, default COB")

    parser_mito.add_argument(
        '--reference', required=False,
        help="Run NOVOplasty in reference mode")

    parser_mito.add_argument(
        '-w', '--workdir', '--tmpdir',
        type=str, dest='workdir',
        required=False,
        help="Temporary directory to store datafiles and processes in")

    parser_mito.add_argument(
        '--pipe', action='store_true',
        help="AAFTF is running in pipeline mode")

    ##########
    # filter
    ##########
    # arguments
    # -i / --indir:  input dir
    # -o / --out: in fastq reads and output prefix
    # -a / --screen_accessions: screening sequence GenBank accessions
    # -u / --screen_urls: screening sequence URLs (fasta format)
    # --debug: print debug messages and do no remove contamdb BAM file
    # read info, either paired data are required or singleton
    # --left: left or forward reads
    # --right: right or reverse reads
    # --single: single unpaired reads
    # or value from --prefix
    # --aligner: bbduk bwa, bowtie2, minimap for read alignment to contamdb

    parser_filter = subparsers.add_parser(
        'filter',
        aliases=['filter_reads', 'read_filter'],
        description="Filter reads which match " +
        "contaminant databases such as phiX",
        help='Filter contaminanting reads')

    parser_filter.add_argument('-w', '--workdir', type=str,
                               help="temp directory")

    parser_filter.add_argument('-c', '--cpus', type=int,
                               metavar="cpus", required=False,
                               default=1,
                               help="Number of CPUs/threads to use.")

    parser_filter.add_argument('-o', '--out',
                               dest='basename', type=str,
                               required=False,
                               help="Output basename")

    parser_filter.add_argument(
        '-v', '--debug', action='store_true',
        help="Provide debugging messages and do not remove contamdb " +
        "matching BAM file")

    parser_filter.add_argument(
        '-a', '--screen_accessions', type=str,
        nargs="*",
        help="Genbank accession number(s) to screen out from initial reads.")

    parser_filter.add_argument(
        '-u', '--screen_urls', type=str,
        nargs="*",
        help="URLs to download and screen out initial reads.")

    parser_filter.add_argument(
        '-s', '--screen_local', type=str,
        nargs="+",
        help="Local FASTA file(s) to use contamination screen")

    parser_filter.add_argument(
        '-l', '--left', required=True,
        help="Left (Forward) reads")

    parser_filter.add_argument('-r', '--right', required=False,
                               help="Right (Reverse) reads")

    parser_filter.add_argument(
        '--AAFTF_DB', type=str,
        required=False,
        help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_filter.add_argument(
        '--aligner', default='bbduk',
        choices=['bbduk', 'bowtie2', 'bwa', 'minimap2'],
        help='Aligner to use to map reads to contamination database')

    parser_filter.add_argument('-m', '--memory', type=int,
                               dest='memory', required=False,
                               help="Max Memory (in GB)")

    parser_filter.add_argument('--pipe', action='store_true',
                               help="AAFTF is running in pipeline mode")

    ##########
    # assemble
    ##########
    # arguments
    # -i / --indir:  input folder
    # -o / --outdir: output folder
    # -p / --prefix: input/outfile prefix
    # --paired or --unpaired
    # --spades
    # --merged for merged reads
    # --tmpdir: tempdir for spades

    parser_asm = subparsers.add_parser(
        'assemble',
        aliases=['asm', 'spades'],
        description="Run assembler on cleaned reads",
        help='Assemble reads')

    parser_asm.add_argument(
        '--method', type=str,
        choices=['spades', 'dipspades', 'megahit', 'masurca', 'nextdenovo'],
        required=False, default="spades",
        help="Assembly method: spades, dipspades, megahit, nextdenovo, masurca")

    parser_asm.add_argument(
        '-o', '--out', type=str,
        required=True,  # think about sensible replacement in future
        help="Output assembly FASTA")

    parser_asm.add_argument('-w', '--workdir', type=str,
                            dest='workdir',
                            help="assembly output directory")

    parser_asm.add_argument('-c', '--cpus', type=int,
                            metavar="cpus",
                            required=False,
                            default=1,
                            help="Number of CPUs/threads to use.")

    parser_asm.add_argument(
        '-m', '--memory', type=str,
        dest='memory',
        required=False,
        default='32',
        help="Memory (in GB) setting for SPAdes. Default is 32gb")

    parser_asm.add_argument('-l', '--left', required=False,
                            help="Left (Forward) reads")

    parser_asm.add_argument('-r', '--right', required=False,
                            help="Right (Reverse) reads")
    parser_asm.add_argument('-lr', '--longreads', required=False,
                            help="Long Read fastq (pacbio or ONT)")
    parser_asm.add_argument('--merged', required=False,
                            help="Merged reads from flash or fastp")

    parser_asm.add_argument('--careful', required=False,
                            action='store_true', default=True,
                            help="Run --careful mode in spades (Default)")

    parser_asm.add_argument('--no-careful', required=False,
                            action='store_false',
                            default=True,
                            dest='careful')

    parser_asm.add_argument(
        '--isolate', required=False,
        action='store_true',
        help="Run --isolate mode not --careful mode in spades")

    parser_asm.add_argument('-v', '--debug',
                            action='store_true',
                            help="Print Spades stdout to terminal")

    parser_asm.add_argument('--tmpdir', type=str,
                            required=False,
                            help="Assembler temporary dir")

    parser_asm.add_argument('--assembler_args',
                            action='append',
                            required=False,
                            help="Additional SPAdes/Megahit arguments")

    parser_asm.add_argument('--haplocontigs',
                            dest='haplocontigs',
                            default=False,
                            action='store_true',
                            help="For dipSPAdes take the haplocontigs file")

    parser_asm.add_argument('--pipe',
                            action='store_true',
                            help="AAFTF is running in pipeline mode")

    ##########
    # vecscreen
    ##########
    # arguments
    # -i / --input:  input assembly file
    # -o / --outfile: output cleaned assembly
    # --prefix: Prefix for output / temp files
    # --tmpdir
    # --pid / percent_id

    parser_vecscreen = subparsers.add_parser(
        'vecscreen',
        aliases=['vectorscreen', 'vector_blast'],
        description="Screen contigs for vector and common contaminantion",
        help='BLASTN Vector and Contaminant Screening of contigs')

    parser_vecscreen.add_argument('-c', '--cpus',
                                  type=int,
                                  metavar="cpus",
                                  default=1,
                                  help="Number of CPUs/threads to use.")

    parser_vecscreen.add_argument(
        '-i', '--input', '--infile',
        type=str,
        required=True,
        dest='infile',
        help="Input contigs or scaffold assembly")

    parser_vecscreen.add_argument(
        '-o', '--outfile',
        type=str,
        required=True,
        help="Output vector screened and cleaned assembly")

    parser_vecscreen.add_argument(
        '-pid', '--percent_id', type=int,
        required=False,
        help="Percent Identity cutoff for vecscreen adaptor matches")

    parser_vecscreen.add_argument(
        '--prefix', type=str,
        required=False, help="Prefix for tempfiles")

    parser_vecscreen.add_argument(
        '-w', '--workdir', '--tmpdir',
        type=str,
        help="Working directory to store datafiles and processes in")

    parser_vecscreen.add_argument(
        '--AAFTF_DB', type=str,
        required=False,
        help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_vecscreen.add_argument(
        '-s', '--stringency',
        default='high', choices=['high', 'low'],
        help="Stringency to filter VecScreen hits")

    parser_vecscreen.add_argument(
        '-v', '--debug', action='store_true', dest='debug',
        help="Provide debugging messages")

    parser_vecscreen.add_argument(
        '--pipe', action='store_true',
        help="AAFTF is running in pipeline mode")

    ##########
    # fcs_screen
    ##########
    # arguments
    # -i / --input:  input assembly file
    # -o / --outfile: output cleaned assembly
    # --prefix: Prefix for output / temp files
    # --euk - expect eukaryotic screening
    # --prol - expect prokaryote screening

    parser_fcs_screen = subparsers.add_parser(
        'fcs_screen',
        aliases=['ncbi_fcs', 'ncbi_fcs-screen'],
        description="Screen with NCBI fcs tool contigs for vector and common contaminantion",
        help='NCBI Foreign Contaminant Screening for Vector sequences in contigs')

    parser_fcs_screen.add_argument(
        '-i', '--input', '--infile',
        type=str,
        required=True,
        dest='infile',
        help="Input contigs or scaffold assembly")

    parser_fcs_screen.add_argument(
        '-o', '--outfile',
        type=str,
        required=True,
        help="Output vector screened and cleaned assembly")

    parser_fcs_screen.add_argument(
        '--prefix', type=str,
        required=False, help="Prefix for tempfiles")

    parser_fcs_screen.add_argument(
        '--container_engine', type=str, default='singularity',
        help="Container engine (singular or docker)")

    parser_fcs_screen.add_argument(
        '--image', type=str, required=False,
        help="Container file (or will download and look in AAFTF_DB)")

    parser_fcs_screen.add_argument(
        '-w', '--workdir', '--tmpdir',
        type=str,
        help="Working directory to store datafiles and processes in")

    parser_fcs_screen.add_argument(
        '--AAFTF_DB', type=str,
        required=False,
        help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_fcs_screen.add_argument(
        '--prok', action='store_true',
        help="Run in Prokaryote matching mode")

    parser_fcs_screen.add_argument(
        '--euk', action='store_true',
        help="Run in Eukaryote matching mode (Default)")

    parser_fcs_screen.add_argument(
        '--fcs_script',  type=str, required=False,
        help="location of the run_fcsadaptor.sh script (or will download automatically)")

    parser_fcs_screen.add_argument(
        '-v', '--debug', action='store_true', dest='debug',
        help="Provide debugging messages")

    parser_fcs_screen.add_argument(
        '--pipe', action='store_true',
        help="AAFTF is running in pipeline mode")

    ##########
    # fcs_gx_purge
    ##########
    # arguments
    # -a / --assembly: input assembly file
    # -o / --out: output cleaned assembly file
    # -p / --prefix: datafile prefix and temp/output file prefix
    # -c / --cpus: number of cpus
    # -d / --db: path to db (this is best if is memory mapped location)
    # -t / --taxid: ncbi taxonid (can be any node in taxonomy)

    parser_fcsgx = subparsers.add_parser(
        'fcs_gx_purge',
        aliases=['ncbi_fcs-gx', 'ncbi_fcs_gx', 'gx'],
        description="Purge contigs based on fcs_gx results",
        help='Purge contigs based on contamination search with fcs_gx')

    parser_fcsgx.add_argument(
        '-i', '--input', type=str,
        required=True,
        help="Input contigs or scaffold assembly")

    parser_fcsgx.add_argument(
        '-o', '--outfile', type=str,
        required=True,  # think about sensible replacement in future
        help="Output fcs_gx cleaned assembly")

    parser_fcsgx.add_argument(
        '--prefix', type=str,
        required=False, help="Prefix for tempfiles")

    parser_fcsgx.add_argument(
        '-t', '--taxid',
        required=False, default=4890,
        help="NCBI Taxonomy ID for contamaination matches, i.e. Ascomycota")

    parser_fcsgx.add_argument(
        '-d', '--db', required=False, default='/my_tmpfs/gxdb/all',
        help="gxdb database path")

    parser_fcsgx.add_argument('-c', '--cpus', type=int,
                              metavar="cpus", default=1,
                              help="Number of CPUs/threads to use.")

    parser_fcsgx.add_argument(
        '-w', '--workdir', '--tmpdir',
        type=str, dest='workdir',
        required=False,
        help="Temporary directory to store datafiles and processes in")

    parser_fcsgx.add_argument('-v', '--debug',
                              action='store_true',
                              dest='debug',
                              help="Provide debugging messages")

    parser_fcsgx.add_argument(
        '--AAFTF_DB', type=str,
        required=False,
        help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_fcsgx.add_argument(
        '--pipe', action='store_true',
        help="AAFTF is running in pipeline mode")

    ##########
    # sourpurge
    ##########
    # arguments
    # -i / --input: input assembly file
    # -o / --out: output cleaned assembly file
    # -p / --prefix: datafile prefix and temp/output file prefix
    # -l / --left: left read of pair
    # -r / --right: right read of pair
    # -k / --kmer: kmer size when running sourmash (needs to match db), default '31'
    # --sourdb_type gtdb gtdbrep gbk
    # --tmpdir
    # --phylum: phylum to keep
    # -mc / --mincovpct:

    parser_sour = subparsers.add_parser(
        'sourpurge',
        aliases=['purge'],
        description="Purge contigs based on sourmash results",
        help='Purge contigs based on sourmash results')

    parser_sour.add_argument(
        '-i', '--input', type=str,
        required=True,
        help="Input contigs or scaffold assembly")

    parser_sour.add_argument(
        '-o', '--outfile', type=str,
        required=True,  # think about sensible replacement in future
        help="Output sourmash cleaned assembly")

    parser_sour.add_argument(
        '--prefix', type=str,
        required=False, help="Prefix for tempfiles")

    parser_sour.add_argument('-l', '--left', required=False,
                             help="Left (Forward) reads")

    parser_sour.add_argument('-r', '--right', required=False,
                             help="Right (Reverse) reads")

    parser_sour.add_argument(
        '-p', '--phylum',
        required=True,
        nargs="+",
        help="Phylum or Phyla to keep matches, i.e. Ascomycota")

    parser_sour.add_argument(
        '--sourdb', required=False,
        help="SourMash LCA taxonomy database (defaults to k-31)")

    parser_sour.add_argument(
        '-k', '--kmer', required=False, default='31',
        help="SourMash LCA kmersize when taxonomy database was built")

    parser_sour.add_argument(
        '-mc', '--mincovpct',
        default=5, type=int,
        help="Minimum percent of N50 coverage to remove")

    parser_sour.add_argument('-c', '--cpus', type=int,
                             metavar="cpus", default=1,
                             help="Number of CPUs/threads to use.")

    parser_sour.add_argument(
        '-w', '--workdir', '--tmpdir',
        type=str, dest='workdir',
        required=False,
        help="Temporary directory to store datafiles and processes in")

    parser_sour.add_argument('-v', '--debug',
                             action='store_true',
                             dest='debug',
                             help="Provide debugging messages")

    parser_sour.add_argument(
        '--sourdb_type', default="gbk",
        required=False, choices=['gbk', 'gtdbrep', 'gtdb'],
        help="Which sourpurge database to use.")

    parser_sour.add_argument(
        '--AAFTF_DB', type=str,
        required=False,
        help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_sour.add_argument(
        '--just-show-taxonomy',
        dest='taxonomy',
        action='store_true',
        help="Show taxonomy information and exit")

    parser_sour.add_argument(
        '--pipe', action='store_true',
        help="AAFTF is running in pipeline mode")

    ##########
    # rmdup
    ##########

    # -i / --input
    # -o / --out
    # --tmpdir
    # --percent_id
    # ---mincovpct
    # -ml / --minlen
    # --exhaustive
    # --debug

    parser_rmdup = subparsers.add_parser(
        'rmdup',
        aliases=['dedup'],
        description="Remove duplicate contigs",
        help='Remove duplicate contigs')
    parser_rmdup.add_argument(
        '-i', '--input', type=str,
        required=True,
        help="Input Assembly fasta file(contigs or scaffolds)")

    parser_rmdup.add_argument(
        '-o', '--out', type=str,
        required=True,
        help="Output new version of assembly with " +
        "duplicated contigs/scaffolds removed")

    parser_rmdup.add_argument('-c', '--cpus', type=int,
                              metavar="cpus", required=False,
                              default=1,
                              help="Number of CPUs/threads to use.")

    parser_rmdup.add_argument(
        '-w', '--workdir',
        '--tmpdir',
        dest='workdir',
        type=str,
        required=False,
        help="Temporary directory to store datafiles and processes in")

    parser_rmdup.add_argument(
        '-pid', '--percent_id', type=int,
        dest='percent_id',
        required=False,
        default=95,
        help="Percent Identity used in matching contigs for redundancy")

    parser_rmdup.add_argument(
        '-pcov', '--percent_cov',
        type=int,
        dest='percent_cov',
        required=False,
        default=95,
        help="Coverage of contig used to decide if it is redundant")

    parser_rmdup.add_argument(
        '-ml', '--minlen', type=int,
        required=False,
        default=500,
        help="Minimum contig length to keep, shorter ones are dropped")

    parser_rmdup.add_argument(
        '--exhaustive', action='store_true',
        help="Compute overlaps for every contig, " +
        "otherwise only process contigs for L75 and below")

    parser_rmdup.add_argument(
        '--debug',
        action='store_true',
        help='Run rmdup in debugging mode for more output')

    parser_rmdup.add_argument('--pipe',
                              action='store_true',
                              help="AAFTF is running in pipeline mode")

    ##########
    # polish
    ##########
    # arguments
    # -i / --in: input assembly file
    # -o / --out: output cleaned assembly
    # --left: left fastq read file (ILLUMINA)
    # --right: right fastq read file (ILLUMINA)
    # --longreads: long reads (ONT/PACBIO) fastq
    # -rp / --reads-prefix: input/outfile reads prefix
    # --method pilon
    # --memory: default 4
    # --iterations: default 5
    # --tmpdir
    # --debug
    # --diploid - inddicate this is a diploid organism

    parser_polish = subparsers.add_parser(
        'polish',
        aliases=['pilon', 'polca'],
        description="Polish contig sequences with Pilon, POLCA, NextPolish",
        help='Polish contig sequences with short reads')

    parser_polish.add_argument('-o', '--out', '--outfile',
                               type=str,
                               dest='outfile', required=False,
                               help="Output a Polished assembly")

    parser_polish.add_argument('-i', '--infile', '--input',
                               type=str, dest='infile',
                               required=True,
                               help="Input contigs or scaffold assembly")

    parser_polish.add_argument('-c', '--cpus', type=int,
                               metavar="cpus", default=1,
                               help="Number of CPUs/threads to use.")

    parser_polish.add_argument('-m', '--memory', type=int,
                               default=16,
                               dest='memory', required=False,
                               help="Max Memory (in GB) (default is 16gb)")

    parser_polish.add_argument('-v', '--debug', action='store_true',
                               help="Provide debugging messages")

    parser_polish.add_argument('-it', '--iterations', type=int,
                               default=5,
                               help="Number of Polishing iterations to run (default is 5)")
    parser_polish.add_argument('--method', type=str,
                               choices=['pilon', 'polca', 'nextpolish', 'racon'],
                               required=False, default="pilon",
                               help="Polishing method: pilon, polca, nextpolish, racon")

    parser_polish.add_argument('--polca', type=str,
                               default="polca.sh",
                               help='polca exe path - provide full path to deal with samtools mismatch in masurca')

    parser_polish.add_argument(
        '-l', '--left', type=str,
        required=False,
        help='The name of the left/forward Illumina reads of paired-end ' +
        'FASTQ formatted reads.')

    parser_polish.add_argument(
        '-r', '--right', type=str,
        required=False,
        help='The name of the right/reverse Illumina reads of paired-end ' +
        'FASTQ formatted reads.')

    parser_polish.add_argument('-lr', '--longreads',
                               required=False,
                               help="Long Read FASTQ (PacBio or ONT)")
    parser_polish.add_argument(
        '-w', '--workdir', '--tmpdir',
        type=str, dest='workdir',
        required=False,
        help="Temporary directory to store datafiles for processes")

    parser_polish.add_argument(
        '--prefix', type=str,
        required=False, help="Prefix for readfiles")

    parser_polish.add_argument(
            '--diploid', action='store_true',
            help="Run pilon in diploid mode - affects heterozygous SNP calling")

    parser_polish.add_argument(
            '--ploidy', default=1, type=int,
            help="Run nextpolish in specific ploidy mode - affects heterozygous SNP calling, default is 1")

    parser_polish.add_argument('--pipe', action='store_true',
                               help="AAFTF is running in pipeline mode")

    ##########
    # sort/rename FASTA headers
    ##########
    # arguments
    # -i / --input: input assembly file
    # -o / --out: output assembly file
    # -n / --name: base name to use default=scaffolds_

    parser_sort = subparsers.add_parser(
        'sort',
        description="Sort contigs by length and rename FASTA headers",
        help='Sort contigs by length and rename FASTA headers')

    parser_sort.add_argument('-i', '--input', '--infile',
                             required=True, dest='input',
                             help='Input genome assembly FASTA')

    parser_sort.add_argument('-o', '--out', '--output',
                             required=True, dest='out',
                             help='Output genome assembly FASTA')

    parser_sort.add_argument(
        '-ml', '--minlen', type=int,
        required=False, default=0,
        help="Minimum contig length to keep, shorter ones are dropped")

    parser_sort.add_argument(
        '-n', '--name', '--basename',
        default='scaffold',
        dest='name',
        help='Basename to rename FASTA headers')

    ##########
    # assess completeness
    ##########
    # arguments
    # -i / --input: input assembly file
    # -r / --report: report file (otherwise stdout)
    # -t / --telomere_monomer: telomere repeat monomer pattern (default=[TAA[C]+])
    # -n / --telomere_n_repeats: Minimum number of repeats of telomere monomer (default=2)

    parser_assess = subparsers.add_parser(
        'assess',
        aliases=['stats'],
        description="Assess completeness of genome assembly",
        help='Assess completeness of genome assembly')

    parser_assess.add_argument(
        '-i', '--input', '--infile',
        required=True,
        help='Input genome assembly to test completeness and ' +
        'provide summary statistics')

    parser_assess.add_argument(
        '-r', '--report', type=str,
        help='Filename to save report information otherwise ' +
        'will print to stdout')

    parser_assess.add_argument(
        '-t', '--telomere_monomer', type=str,
        help='Telomere repeat monomer to search for. default(TAA[C]+)',
        default='TAA[C]+')

    parser_assess.add_argument(
        '-n', '--telomere_n_repeat', type=int, default=2,
        help='Telomere minimum number of monomer repeats. (default 2)')

    ##########
    # pipeline run it all
    ##########
    # arguments
    # -i / --input: input assembly file
    # -r / --report: report file (otherwise stdout)
    # --tmpdir

    parser_pipeline = subparsers.add_parser(
        'pipeline',
        description="Run entire AAFTF pipeline automagically",
        help='Run AAFTF pipeline')

    parser_pipeline.add_argument('--tmpdir', type=str,
                                 required=False,
                                 help="Assembler temporary dir")
    parser_pipeline.add_argument('--assembler_args', action='append',
                                 required=False,
                                 help="Additional SPAdes/Megahit arguments")
    parser_pipeline.add_argument(
        '--method', type=str,
        required=False, default="spades",
        help="Assembly method: spades, dipspades, megahit")

    parser_pipeline.add_argument(
        '-l', '--left', type=str,
        required=True,
        help='left/forward reads of paired-end FASTQ or ' +
        'single-end FASTQ.')

    parser_pipeline.add_argument(
        '-r', '--right', type=str,
        required=False,
        help='right/reverse reads of paired-end FASTQ.')

    parser_pipeline.add_argument(
        '-o', '--out', type=str,
        required=True, dest='basename',
        help="Output basename, default to base name of --left reads")

    parser_pipeline.add_argument('-c', '--cpus', type=int, metavar="cpus",
                                 required=False, default=1,
                                 help="Number of CPUs/threads to use.")

    parser_pipeline.add_argument(
        '-m', '--memory', type=str,
        dest='memory', required=False,
        help="Memory (in GB) setting for SPAdes. Default is Auto")

    parser_pipeline.add_argument(
        '-ml', '--minlen', type=int,
        default=75,
        required=False,
        help="Minimum read length after trimming, default: 75")

    parser_pipeline.add_argument(
        '-a', '--screen_accessions',
        type=str,
        nargs="*",
        help="Genbank accession number(s) to screen out from initial reads.")

    parser_pipeline.add_argument(
        '-u', '--screen_urls', type=str,
        nargs="*",
        help="URLs to download and screen out initial reads.")

    parser_pipeline.add_argument(
        '-it', '--iterations',
        type=int, default=5,
        help="Number of Pilon Polishing iterations to run")

    parser_pipeline.add_argument('-mc', '--mincontiglen', type=int,
                                 default=500,
                                 required=False,
                                 help="Minimum length of contigs to keep")

    parser_pipeline.add_argument(
        '--AAFTF_DB', type=str,
        required=False,
        help='Path to AAFTF resources, defaults to $AAFTF_DB')

    parser_pipeline.add_argument('-w', '--workdir', type=str,
                                 help="temp directory")

    parser_pipeline.add_argument('-v', '--debug', action='store_true',
                                 help="Provide debugging messages")

    parser_pipeline.add_argument(
        '-p', '--phylum', required=True,
        nargs="+",
        help="Phylum or Phyla to keep matches, i.e. Ascomycota")

    parser_pipeline.add_argument(
        '--sourdb', required=False,
        help="SourMash LCA k-31 taxonomy database")

    parser_pipeline.add_argument(
        '--mincovpct', default=5, type=int,
        help="Minimum percent of N50 coverage to remove")

    # set defaults
    parser.set_defaults(func=run_subtool)

    # process the arguments now
    # if no args then print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    try:
        status(f'Running AAFTF v{myversion}')
        args.func(parser, args)
    except OSError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


if __name__ == "__main__":
    main()
