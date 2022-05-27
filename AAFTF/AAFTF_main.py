#!/usr/bin/env python3

# note structure of code taken from poretools https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import os.path
import sys
import argparse

# AAFTF imports
from AAFTF.__version__ import __version__
myversion = __version__
from AAFTF.utility import status

def run_subtool(parser, args):
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
    elif args.command == 'sourpurge':
        import AAFTF.sourpurge as submodule
    elif args.command == 'rmdup':
        import AAFTF.rmdup as submodule
    elif args.command == 'pilon':
        import AAFTF.pilon as submodule
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

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                            action="store_true",
                            dest="quiet")
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():

    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='AAFTF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")
    parser.add_argument("-v", "--version", help="Installed AAFTF version",
                        action="version",
                        version="%(prog)s " + str(myversion))

    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)


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

    # read info, either paired data are required or singleton
    # --left: left or forward reads
    # --right: right or reverse reads
    #
    # --merge: merge reads in fastp

    parser_trim = subparsers.add_parser('trim',
       description="This comamnd trims reads in FASTQ format to remove low quality reads and trim adaptor sequences",
       help='Trim FASTQ input reads')

    parser_trim.add_argument('-o','--out',type=str,
                             required=False, dest='basename',
                             help="Output basename, default to base name of --left reads")

    parser_trim.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                              help="Number of CPUs/threads to use.")

    parser_trim.add_argument('-ml','--minlen',type=int,
                             default=75,
                             required=False,
                             help="Minimum read length after trimming, default: 75")

    parser_trim.add_argument('-l', '--left',type=str,
                              required=True,
            help='left/forward reads of paired-end FASTQ or single-end FASTQ.')

    parser_trim.add_argument('-r', '--right',type=str,
                              required=False,
            help='right/reverse reads of paired-end FASTQ.')

    parser_trim.add_argument('-v','--debug',action='store_true',
                             help="Provide debugging messages")

    parser_trim.add_argument('--pipe',action='store_true',
                             help="AAFTF is running in pipeline mode")

    parser_trim.add_argument('--method', default='bbduk',
                            choices=['bbduk', 'trimmomatic', 'fastp'],
                            help='Program to use for adapter trimming')

    parser_trim.add_argument('-m','--memory',type=int,
                            dest='memory',required=False,
                            help="Max Memory (in GB)")

    parser_trim.add_argument('--merge',action='store_true',
                            help="Merge paired end reads when running fastp")

    tool_group = parser_trim.add_mutually_exclusive_group(required=False)

    tool_group.add_argument('--trimmomatic','--jar', metavar='trimmomatic_jar',
                            type=str,required=False,
                            help='Trimmomatic JAR path')
    trimmomatic_group = parser_trim.add_argument_group(title='Trimmomatic options',
                                              description="Trimmomatic trimming options")

    trimmomatic_group.add_argument('--trimmomatic_adaptors',
                                   default="TruSeq3-PE.fa",
                                   help="Trimmomatic adaptor file, default: TruSeq3-PE.fa")

    trimmomatic_group.add_argument('--trimmomatic_clip',
                                   default="ILLUMINACLIP:%s:2:30:10",
                                   help="Trimmomatic clipping, default: ILLUMINACLIP:TruSeq3-PE.fa:2:30:10")

    trimmomatic_group.add_argument('--trimmomatic_leadingwindow',
                                   default="3",type=int,
                                   help="Trimmomatic window processing arguments, default: LEADING:3")

    trimmomatic_group.add_argument('--trimmomatic_trailingwindow',
                                   default="3",type=int,
                                   help="Trimmomatic window processing arguments, default: TRAILING:3")

    trimmomatic_group.add_argument('--trimmomatic_slidingwindow',
                                   default="4:15",type=str,
                                   help="Trimmomatic window processing arguments, default: SLIDINGWINDOW:4:15")
    trimmomatic_group.add_argument('--trimmomatic_quality',
                                   default="phred33",
                                   help="Trimmomatic quality encoding -phred33 or phred64")



    ##########
    # mito-asm assembly mitochondrial genome
    ##########

    parser_mito = subparsers.add_parser('mito',
                                        description="De novo assembly of mitochondrial genome using NOVOplasty, takes PE Illumina adapter trimmed data.",
                                        help='De novo assembly of mitochondrial genome')
    parser_mito.add_argument('-l', '--left',required=True,
                             help="Left (Forward) reads")

    parser_mito.add_argument('-r', '--right',required=True,
                             help="Right (Reverse) reads")

    parser_mito.add_argument('-o','--out',type=str,
                            required=True,
                            help="Output FASTA file for mitochondrial genome")

    parser_mito.add_argument('--minlen',default=10000,type=int,
                             help="Minimum expected genome size")

    parser_mito.add_argument('--maxlen',default=100000,type=int,
                             help="Maximum expected genome size")

    parser_mito.add_argument('-s','--seed',required=False,
                             help="Seed sequence, ie related mitochondrial genome, Default: A. nidulans")

    parser_mito.add_argument('--starting',required=False,
                             help="FASTA file of start sequence, rotate genome to, default COB")

    parser_mito.add_argument('--reference',required=False,
                             help="Run NOVOplasty in reference mode")

    parser_mito.add_argument('-w', '--workdir', '--tmpdir',
                            type=str, dest='workdir',
                            required=False,
                            help="Temporary directory to store datafiles and processes in")

    parser_mito.add_argument('--pipe',action='store_true',
                            help="AAFTF is running in pipeline mode")


    ##########
    # filter
    ##########
    # arguments
    # -i / --indir:  input dir
    # -p / --prefix: in fastq reads and output prefix
    # -a / --screen_accessions: screening sequence GenBank accessions
    # -u / --screen_urls: screening sequence URLs (fasta format)
    # --debug: print debug messages and do no remove contamdb BAM file
    # read info, either paired data are required or singleton
    # --left: left or forward reads
    # --right: right or reverse reads
    # --single: single unpaired reads
    # or value from --prefix
    # --aligner: bbduk bwa, bowtie2, minimap for read alignment to contamdb

    parser_filter = subparsers.add_parser('filter',
        description="Filter reads which match contaminant databases such as phiX",
    help='Filter contaminanting reads')

    parser_filter.add_argument('-w', '--workdir',type=str,
                        help="temp directory")

    parser_filter.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                        help="Number of CPUs/threads to use.")

    parser_filter.add_argument('-o','--out',dest='basename', type=str,
                        required=False,
                        help="Output basename")

    parser_filter.add_argument('-v','--debug',action='store_true',
                             help="Provide debugging messages and do not remove contamdb matching BAM")

    parser_filter.add_argument('-a','--screen_accessions',type = str,
                               nargs="*",
                               help="Genbank accession number(s) to screen out from initial reads.")

    parser_filter.add_argument('-u','--screen_urls',type = str,
                               nargs="*",
                               help="URLs to download and screen out initial reads.")

    parser_filter.add_argument('-s','--screen_local',type = str,
                               nargs="+",
                               help="Local FASTA file(s) to use contamination screen")

    parser_filter.add_argument('-l', '--left',required=True,
                             help="Left (Forward) reads")

    parser_filter.add_argument('-r', '--right',required=False,
                             help="Right (Reverse) reads")

    parser_filter.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_filter.add_argument('--aligner', default='bbduk',
                               choices=['bbduk', 'bowtie2', 'bwa', 'minimap2'],
                               help='Aligner to use to map reads to contamination database')

    parser_filter.add_argument('-m','--memory',type=int,
                            dest='memory',required=False,
                            help="Max Memory (in GB)")

    parser_filter.add_argument('--pipe',action='store_true',
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

    parser_asm = subparsers.add_parser('assemble',
                                       description="Run assembler on cleaned reads",
                                       help='Assemble reads')

    parser_asm.add_argument('--method',type=str,
                             required=False, default="spades",
                             help="Assembly method: spades, dipspades, megahit")

    parser_asm.add_argument('-o','--out',type=str,
                             required=True, # think about sensible replacement in future
                             help="Output assembly FASTA")

    parser_asm.add_argument('-w', '--workdir',type=str,
                        dest='workdir',
                        help="assembly output directory")

    parser_asm.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                        help="Number of CPUs/threads to use.")

    parser_asm.add_argument('-m','--memory',type=str,
                            dest='memory',required=False,default='32',
                            help="Memory (in GB) setting for SPAdes. Default is 32")

    parser_asm.add_argument('-l', '--left',required=False,
                             help="Left (Forward) reads")

    parser_asm.add_argument('-r', '--right',required=False,
                             help="Right (Reverse) reads")

    parser_asm.add_argument('--merged',required=False,
                             help="Merged reads from flash or fastp")

    parser_asm.add_argument('-v','--debug',action='store_true',
                             help="Print Spades stdout to terminal")

    parser_asm.add_argument('--tmpdir',type=str,required=False,help="Assembler temporary dir")
    parser_asm.add_argument('--assembler_args',action='append',required=False,help="Additional SPAdes/Megahit arguments")
    parser_asm.add_argument('--haplocontigs',dest='haplocontigs',default=False, action='store_true',help="For dipSPAdes take the haplocontigs file")

    parser_asm.add_argument('--pipe',action='store_true',
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

    parser_vecscreen = subparsers.add_parser('vecscreen',
                                             description="Screen contigs for vector and common contaminantion",
                                             help='Vector and Contaminant Screening of assembled contigs')

    parser_vecscreen.add_argument('-c','--cpus',type=int,metavar="cpus",default=1,
                                  help="Number of CPUs/threads to use.")

    parser_vecscreen.add_argument('-i','--input','--infile',type=str,
                                  required=True, dest='infile',
                                  help="Input contigs or scaffold assembly")

    parser_vecscreen.add_argument('-o','--outfile',type=str,
                                  required=True,
                                  help="Output vector screened and cleaned assembly")

    parser_vecscreen.add_argument('-pid','--percent_id',type=int,
                                  required=False,
                                  help="Percent Identity cutoff for vecscreen adaptor matches")

    parser_vecscreen.add_argument('-w', '--workdir', '--tmpdir',type=str,
                        help="Working directory to store datafiles and processes in")

    parser_vecscreen.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_vecscreen.add_argument('-s', '--stringency', default='high', choices=['high','low'],
                                  help="Stringency to filter VecScreen hits")

    parser_vecscreen.add_argument('-v','--debug', action='store_true', dest='debug',
                             help="Provide debugging messages")

    parser_vecscreen.add_argument('--pipe',action='store_true',
                             help="AAFTF is running in pipeline mode")


    ##########
    # sourpurge
    ##########
    # arguments
    # -a / --assembly: input assembly file
    # -o / --out: output cleaned assembly file
    # -p / --prefix: datafile prefix and temp/output file prefix
    # -i / --indir: directory where sequence reads are located
    # -c / --cpus: number of cpus
    # --tmpdir
    # --phylum: phylum to keep
    parser_sour = subparsers.add_parser('sourpurge',
                                        description="Purge contigs based on sourmash results",
                                        help='Purge contigs based on sourmash results')

    parser_sour.add_argument('-i','--input',type=str,
                             required=True,
                             help="Input contigs or scaffold assembly")

    parser_sour.add_argument('-o','--outfile',type=str,
                             required=True, # think about sensible replacement in future
                             help="Output sourmash cleaned assembly")

    parser_sour.add_argument('-l', '--left',required=False,
                             help="Left (Forward) reads")

    parser_sour.add_argument('-r', '--right',required=False,
                             help="Right (Reverse) reads")

    parser_sour.add_argument('-p', '--phylum',required=True, nargs="+",
                             help="Phylum or Phyla to keep matches, i.e. Ascomycota")

    parser_sour.add_argument('--sourdb',required=False,
                             help="SourMash LCA k-31 taxonomy database")

    parser_sour.add_argument('-mc', '--mincovpct',default=5,type=int,
                             help="Minimum percent of N50 coverage to remove")

    parser_sour.add_argument('-c','--cpus',type=int,metavar="cpus",default=1,
                                  help="Number of CPUs/threads to use.")

    parser_sour.add_argument('-w', '--workdir', '--tmpdir',type=str, dest='workdir',
                        required=False,
                        help="Temporary directory to store datafiles and processes in")

    parser_sour.add_argument('-v','--debug', action='store_true', dest='debug',
                             help="Provide debugging messages")

    parser_sour.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_sour.add_argument('--just-show-taxonomy',dest='taxonomy', action='store_true',
                               help="Show taxonomy information and exit")

    parser_sour.add_argument('--pipe',action='store_true',
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

    parser_rmdup = subparsers.add_parser('rmdup',
                                         description="Remove duplicate contigs",
                                         help='Remove duplicate contigs')
    parser_rmdup.add_argument('-i','--input',type=str,
                               required=True,
                               help="Input Assembly fasta file(contigs or scaffolds)")

    parser_rmdup.add_argument('-o','--out',type=str,
                               required=True,
                               help="Output new version of assembly with duplicated contigs/scaffolds removed")

    parser_rmdup.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                        help="Number of CPUs/threads to use.")

    parser_rmdup.add_argument('-w', '--workdir', '--tmpdir', dest='workdir',type=str,
                               required=False,
                               help="Temporary directory to store datafiles and processes in")

    parser_rmdup.add_argument('-pid','--percent_id',type=int, dest='percent_id',
                               required=False,default=95,
                               help="Percent Identity used in matching contigs for redundancy")

    parser_rmdup.add_argument('-pcov','--percent_cov',type=int, dest='percent_cov',
                               required=False,default=95,
                               help="Coverage of contig used to decide if it is redundant")

    parser_rmdup.add_argument('-ml','--minlen',type=int,
                               required=False,default=500,
                               help="Minimum contig length to keep, shorter ones are dropped")

    parser_rmdup.add_argument('--exhaustive',action='store_true',
                               help="Compute overlaps for every contig, otherwise only process contigs for L75 and below")

    parser_rmdup.add_argument('--debug',action='store_true', help='Run rmdup in debugging mode for more output')

    parser_rmdup.add_argument('--pipe',action='store_true',
                             help="AAFTF is running in pipeline mode")

    ##########
    # pilon
    ##########
    # arguments
    # -i / --in: input assembly file
    # -o / --out: output cleaned assembly
    # -rp / --reads-prefix: input/outfile reads prefix
    # --memory: default 4
    # --iterations: default 5
    # --tmpdir
    # --debug

    parser_pilon = subparsers.add_parser('pilon',
                                         description="Polish contig sequences with Pilon",
                                         help='Polish contig sequences with Pilon')

    parser_pilon.add_argument('-o','--out','--outfile', type=str, dest='outfile',
                             required=True,
                             help="Output Pilon polished assembly")

    parser_pilon.add_argument('-i','--infile','--input', type=str, dest='infile',
                              required=True,
                              help="Input contigs or scaffold assembly")

    parser_pilon.add_argument('-c','--cpus',type=int,metavar="cpus",default=1,
                                  help="Number of CPUs/threads to use.")

    parser_pilon.add_argument('-m','--memory',type=int,default=16,
                            dest='memory',required=False,
                            help="Max Memory (in GB) (default is 16gb)")

    parser_pilon.add_argument('-v','--debug',action='store_true',
                              help="Provide debugging messages")

    parser_pilon.add_argument('-it','--iterations', type=int, default=5,
                              help="Number of Polishing iterations to run (default is 5)")

    parser_pilon.add_argument('-l', '--left',type=str,
                              required=True,
            help='The name of the left/forward reads of paired-end FASTQ formatted reads.')

    parser_pilon.add_argument('-r', '--right',type=str,
                              required=True,
            help='The name of the right/reverse reads of paired-end FASTQ formatted reads.')

    parser_pilon.add_argument('-w', '--workdir', '--tmpdir',
                              type=str, dest='workdir',
                              required=False,
                              help="Temporary directory to store datafiles and processes in")

    parser_pilon.add_argument('--pipe',action='store_true',
                             help="AAFTF is running in pipeline mode")

    ##########
    # sort/rename FASTA headers
    ##########
    # arguments
    # -i / --input: input assembly file
    # -o / --out: output assembly file
    # -n / --name: base name to use default=scaffolds_

    parser_sort = subparsers.add_parser('sort',
                                         description="Sort contigs by length and rename FASTA headers",
                                         help='Sort contigs by length and rename FASTA headers')

    parser_sort.add_argument('-i','--input','--infile',required=True, dest='input',
                               help='Input genome assembly FASTA')

    parser_sort.add_argument('-o','--out','--output',required=True, dest='out',
                               help='Output genome assembly FASTA')

    parser_sort.add_argument('-ml','--minlen',type=int,
                             required=False,default=0,
                             help="Minimum contig length to keep, shorter ones are dropped")

    parser_sort.add_argument('-n','--name','--basename',default='scaffold', dest='name',
                               help='Basename to rename FASTA headers')

    ##########
    # assess completeness
    ##########
    # arguments
    # -i / --input: input assembly file
    # -r / --report: report file (otherwise stdout)
    # --tmpdir

    parser_assess = subparsers.add_parser('assess',
                                          description="Assess completeness of genome assembly",
                                          help='Assess completeness of genome assembly')

    parser_assess.add_argument('-i','--input','--infile',required=True,
                               help='Input genome assembly to test completeness and provide summary statistics')

    parser_assess.add_argument('-r','--report',type=str,
                               help='Filename to save report information otherwise will print to stdout')


    ##########
    # pipeline run it all
    ##########
    # arguments
    # -i / --input: input assembly file
    # -r / --report: report file (otherwise stdout)
    # --tmpdir

    parser_pipeline = subparsers.add_parser('pipeline',
                            description="Run entire AAFTF pipeline automagically",
                            help='Run AAFTF pipeline')

    parser_pipeline.add_argument('--tmpdir',type=str,required=False,help="Assembler temporary dir")
    parser_pipeline.add_argument('--assembler_args',action='append',required=False,help="Additional SPAdes/Megahit arguments")
    parser_pipeline.add_argument('--method',type=str,
                             required=False, default="spades",
                             help="Assembly method: spades, dipspades, megahit")

    parser_pipeline.add_argument('-l', '--left',type=str,
                              required=True,
            help='left/forward reads of paired-end FASTQ or single-end FASTQ.')

    parser_pipeline.add_argument('-r', '--right',type=str,
                              required=False,
            help='right/reverse reads of paired-end FASTQ.')

    parser_pipeline.add_argument('-o','--out',type=str,
                             required=True, dest='basename',
                             help="Output basename, default to base name of --left reads")

    parser_pipeline.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                              help="Number of CPUs/threads to use.")

    parser_pipeline.add_argument('-m','--memory',type=str,
                            dest='memory',required=False,
                            help="Memory (in GB) setting for SPAdes. Default is Auto")

    parser_pipeline.add_argument('-ml','--minlen',type=int,
                             default=75,
                             required=False,
                             help="Minimum read length after trimming, default: 75")

    parser_pipeline.add_argument('-a','--screen_accessions',type = str,
                               nargs="*",
                               help="Genbank accession number(s) to screen out from initial reads.")

    parser_pipeline.add_argument('-u','--screen_urls',type = str,
                               nargs="*",
                               help="URLs to download and screen out initial reads.")

    parser_pipeline.add_argument('-it','--iterations', type=int, default=5,
                              help="Number of Pilon Polishing iterations to run")

    parser_pipeline.add_argument('-mc','--mincontiglen',type=int,
                             default=500,
                             required=False,
                             help="Minimum length of contigs to keep")

    parser_pipeline.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_pipeline.add_argument('-w', '--workdir',type=str,
                        help="temp directory")

    parser_pipeline.add_argument('-v','--debug',action='store_true',
                             help="Provide debugging messages")

    parser_pipeline.add_argument('-p', '--phylum',required=True, nargs="+",
                             help="Phylum or Phyla to keep matches, i.e. Ascomycota")

    parser_pipeline.add_argument('--sourdb',required=False,
                             help="SourMash LCA k-31 taxonomy database")

    parser_pipeline.add_argument('--mincovpct',default=5,type=int,
                             help="Minimum percent of N50 coverage to remove")


    #set defaults
    parser.set_defaults(func=run_subtool)

    ### process args now ###
    # if no args then print help and exit
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    try:
        status('Running AAFTF v{:}'.format(myversion))
        args.func(parser, args)
    except IOError as e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
