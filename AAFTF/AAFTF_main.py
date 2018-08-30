#!/usr/bin/env python3

# note structure of code taken from poretools https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import os.path
import sys
import argparse

import logging
logger = logging.getLogger('AAFTF')

# AAFTF imports
from AAFTF.version import __version__
myversion = __version__

def run_subtool(parser, args):
    if args.command == 'trim':
        import AAFTF.trim as submodule
    elif args.command == 'filter':
        import AAFTF.filter as submodule
    elif args.command == 'assemble':
        import AAFTF.assemble as submodule
    elif args.command == 'vecscreen':
        import AAFTF.vecscreen as submodule
    elif args.command == 'blobpurge':
        import AAFTF.blobpurge as submodule
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
    logging.basicConfig()

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
    # --trimmomatic: arguments are path to JAR or application respectively
    # assume java is PATH already for trimmomatic
    # -o / --outdir: write outdir
    # -p / --prefix: outfile prefix
    # -ml / --minlength: min read length
    
    # read info, either paired data are required or singleton
    # --left: left or forward reads
    # --right: right or reverse reads
    # currently singleton / unpaired reads not supported?
    
    parser_trim = subparsers.add_parser('trim',
       description="This comamnd trims reads in FASTQ format to remove low quality reads and trim adaptor sequences",
       help='Trim FASTQ input reads')
    
    parser_trim.add_argument('-p','--prefix',type=str,
                             required=False,
                             help="Output Prefix, default to base name of --left reads")

    parser_trim.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                              help="Number of CPUs/threads to use.")
    
    parser_trim.add_argument('-o','--outdir', '-w', '--workdir',type=str,
                             default='working_AAFTF', dest='outdir',
                             help="Output directory for trimmed reads")

    parser_trim.add_argument('-ml','--minlength',type=int,
                             default=75,
                             required=False,
                             help="Minimum read length after trimming, default: 75")
    
    parser_trim.add_argument('--left',type=str,
                              required=True,
            help='left/forward reads of paired-end FASTQ or single-end FASTQ.')

    parser_trim.add_argument('--right',type=str,
                              required=False,
            help='right/reverse reads of paired-end FASTQ.')

    parser_trim.add_argument('-v','--debug',action='store_true',
                             help="Provide debugging messages")

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
    # or value from --prefix 
    # --aligner: bwa, bowtie2, minimap for read alignment to contamdb

    parser_filter = subparsers.add_parser('filter',
        description="Filter reads which match contaminant databases such as phiX",
    help='Filter contaminanting reads')

    parser_filter.add_argument('-w', '--workdir',type=str,
                        default="working_AAFTF",
                        help="Working directory to store datafiles and processes in")

    parser_filter.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                        help="Number of CPUs/threads to use.")
    
    parser_filter.add_argument('-p','--prefix',type=str,
                        required=False,
                        help="Input/Output Prefix for fileset")

    parser_filter.add_argument('-v','--debug',action='store_true',
                             help="Provide debugging messages and do not remove contamdb matching BAM")

    parser_filter.add_argument('-a','--screen_accessions',type = str,
                               nargs="*",
                               help="Genbank accession number(s) to screen out from initial reads.")
                               
    parser_filter.add_argument('-u','--screen_urls',type = str,
                               nargs="*",
                               help="URLs to download and screen out initial reads.")

    parser_filter.add_argument('--left',required=False,
                             help="Left (Forward) reads")

    parser_filter.add_argument('--right',required=False,
                             help="Right (Reverse) reads")
    
    parser_filter.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")
    
    parser_filter.add_argument('--aligner', default='bwa', 
                               choices=['bowtie2', 'bwa', 'minimap2'],
                               help='Aligner to use to map reads to contamination database')
    

    ##########
    # assemble
    ##########
    # arguments
    # -i / --indir:  input folder
    # -o / --outdir: output folder
    # -p / --prefix: input/outfile prefix
    # --paired or --unpaired
    # --spades

    parser_asm = subparsers.add_parser('assemble',
                                       description="Run assembler on cleaned reads",
                                       help='Assemble reads')
    
    parser_asm.add_argument('-o','--out',type=str,
                             required=False, # think about sensible replacement in future
                             help="Output spades assembly")

    parser_asm.add_argument('-w', '--workdir', '--tmpdir',type=str,
                        dest='workdir' ,default="working_AAFTF",
                        help="Temporary directory to store datafiles and processes in")

    parser_asm.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                        help="Number of CPUs/threads to use.")
    
    parser_asm.add_argument('-p','--prefix',type=str,
                        required=False,
                        help="Input/Output Prefix for fileset")

    parser_asm.add_argument('-m','--memory',type=str,
                            dest='memory',required=False,default='32',
                            help="Memory (in GB) setting for SPAdes. Default is 32")

    parser_asm.add_argument('--left',required=False,
                             help="Left (Forward) reads")

    parser_asm.add_argument('--right',required=False,
                             help="Right (Reverse) reads")

    parser_asm.add_argument('-v','--debug',action='store_true',
                             help="Print Spades stdout to terminal")

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
                                  required=False,
                                  help="Output vector screened and cleaned assembly (defaults to infile.clean.fasta)")

    parser_vecscreen.add_argument('-p','--prefix',type=str,
                                  required=False,
                                  help="Input/Output Prefix for fileset and tempfiles")

    parser_vecscreen.add_argument('-pid','--percent_id',type=int,
                                  required=False,
                                  help="Percent Identity cutoff for vecscreen adaptor matches")

    parser_vecscreen.add_argument('-w', '--workdir', '--tmpdir',type=str,
                        default="working_AAFTF",
                        help="Working directory to store datafiles and processes in")

    parser_vecscreen.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_vecscreen.add_argument('-s', '--stringency', default='high', choices=['high','low'],
                                  help="Stringency to filter VecScreen hits")


    ##########
    # blobpurge
    ##########
    # arguments
    # -a / --assembly: input assembly file
    # -o / --out: output cleaned assembly file
    # -p / --prefix: sequence reads prefix
    # -i / --indir: directory where sequence reads are located
    # -c / --cpus: number of cpus
    # --tmpdir
    # --phylum: phylum to keep
    parser_blob = subparsers.add_parser('blobpurge',
                                        description="Purge contigs based on BlobPlot results",
                                        help='Purge contigs based on BlobPlot results')

    parser_blob.add_argument('-i','--input',type=str,
                             required=True,
                             help="Input contigs or scaffold assembly")

    parser_blob.add_argument('-o','--out',type=str,
                             required=True, # think about sensible replacement in future
                             help="Output blobplot cleaned assembly")

    parser_blob.add_argument('-p','--prefix',required=False,
                             help="Prefix of the sequence reads files")

    parser_blob.add_argument('--left',required=False,
                             help="Left (Forward) reads")

    parser_blob.add_argument('--right',required=False,
                             help="Right (Reverse) reads")

    parser_blob.add_argument('--phylum',required=True,nargs="+",
                             help="Phylum or Phyla to keep matches from megablast")
    
    parser_blob.add_argument('--blastdb',required=True,
                             help="NCBI nt blast db for classifying contigs/scaffolds by taxa")

    parser_blob.add_argument('-c','--cpus',type=int,metavar="cpus",default=1,
                                  help="Number of CPUs/threads to use.")

    parser_blob.add_argument('-e','--evalue',type=str,default="1e-25",
                             help="Megablast e-value cutoff")

    parser_blob.add_argument('--tmpdir',type=str,
                        required=False,default="working_AAFTF",
                        help="Temporary directory to store datafiles and processes in")

    parser_blob.add_argument('-v','--debug',action='store_true',
                             help="Provide debugging messages")
    # remote or local megablast?

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

    parser_sour.add_argument('-p','--prefix',required=False,
                             help="Prefix of the sequence reads files")

    parser_sour.add_argument('--left',required=False,
                             help="Left (Forward) reads")

    parser_sour.add_argument('--right',required=False,
                             help="Right (Reverse) reads")

    parser_sour.add_argument('--phylum',required=True, nargs="+",
                             help="Phylum or Phyla to keep matches, i.e. Ascomycota")
    
    parser_sour.add_argument('--sourdb',required=False,
                             help="SourMash LCA k-31 taxonomy database")

    parser_sour.add_argument('-m', '--mincovpct',default=5,type=int,
                             help="Minimum percent of N50 coverage to remove")

    parser_sour.add_argument('-c','--cpus',type=int,metavar="cpus",default=1,
                                  help="Number of CPUs/threads to use.")

    parser_sour.add_argument('-w', '--workdir', '--tmpdir',type=str, dest='workdir',
                        required=False,default="working_AAFTF",
                        help="Temporary directory to store datafiles and processes in")

    parser_sour.add_argument('-v','--debug',action='store_true',
                             help="Provide debugging messages")

    parser_sour.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")

    parser_sour.add_argument('--just-show-taxonomy',dest='taxonomy', action='store_true',
                               help="Show taxonomy information and exit")

        
    ##########
    # rmdup
    ##########

    # -i / --input
    # -o / --out
    # --tmpdir
    # --percent_id
    # -c / --coverage
    # -m / --minlen
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

    parser_rmdup.add_argument('-p','--prefix',type=str,
                               required=False,
                               help="Prefix of output file names or temp files")
    parser_rmdup.add_argument('-c','--cpus',type=int,metavar="cpus",required=False,default=1,
                        help="Number of CPUs/threads to use.")
    
    parser_rmdup.add_argument('-w', '--workdir', '--tmpdir', dest='workdir',type=str,
                               required=False,default="working_AAFTF",
                               help="Temporary directory to store datafiles and processes in")

    parser_rmdup.add_argument('-pid','--percent_id',type=int, dest='percent_id', 
                               required=False,default=95,
                               help="Percent Identity used in matching contigs for redundancy")

    parser_rmdup.add_argument('-pcov','--percent_cov',type=int, dest='percent_cov',
                               required=False,default=95,
                               help="Coverage of contig used to decide if it is redundant")

    parser_rmdup.add_argument('-m','--minlen',type=int,
                               required=False,default=500,
                               help="Minimum contig length to keep, shorter ones are dropped")

    parser_rmdup.add_argument('--exhaustive',action='store_true',
                               help="Compute overlaps for every contig, otherwise only process contigs for L50 and below")

    parser_rmdup.add_argument('--debug',action='store_true', help='Run rmdup in debugging mode for more output')

    ##########
    # pilon
    ##########
    # arguments
    # -i / --in: input assembly file
    # -o / --out: output cleaned assembly
    # -rp / --reads-prefix: input/outfile reads prefix
    # --iterations: default 5
    # --tmpdir
    parser_pilon = subparsers.add_parser('pilon',
                                         description="Polish contig sequences with Pilon",
                                         help='Polish contig sequences with Pilon')

    parser_pilon.add_argument('-o','--out','--outfile', type=str, dest='outfile', 
                             required=False,
                             help="Output Pilon polished assembly (defaults to infile.pilon.fasta)")

    parser_pilon.add_argument('-i','--infile','--input', type=str, dest='infile',
                              required=True,
                              help="Input contigs or scaffold assembly")

    parser_pilon.add_argument('-c','--cpus',type=int,metavar="cpus",default=1,
                                  help="Number of CPUs/threads to use.")

    parser_pilon.add_argument('-p','--prefix',required=False,
                              help="Prefix of the read pairs ")

    parser_pilon.add_argument('-it','--iterations', type=int, default=5,
                              help="Number of Polishing iterations to run")

    parser_pilon.add_argument('--left',type=str,
                              required=False,
            help='The name of the left/forward reads of paired-end FASTQ formatted reads.')

    parser_pilon.add_argument('--right',type=str,
                              required=False,
            help='The name of the right/reverse reads of paired-end FASTQ formatted reads.')

    parser_pilon.add_argument('-w', '--workdir', '--tmpdir',
                              type=str, dest='workdir',
                              required=False,default="working_AAFTF",
                              help="Temporary directory to store datafiles and processes in")

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

    parser_sort.add_argument('-n','--name','--basename',default='scaffold', dest='name',
                               help='Basename to rename FASTA headers')
  
    ##########
    # assess completeness
    ##########
    # arguments
    # -i / --input: input assembly file
    # -r / --report: report file (otherwise stdout)
    # --stats: only present the summary stats
    # --method: busco or FGMP (eventually)
    # --clade: for busco - path to odb9 folder 
    # --tmpdir

    parser_assess = subparsers.add_parser('assess',
                                         description="Assess completeness of genome assembly",
                                         help='Assess completeness of genome assembly')

    parser_assess.add_argument('-i','--input','--infile',required=True,
                               help='Input genome assembly to test completeness and provide summary statistics')

    parser_assess.add_argument('-r','--report',type=str,
                               help='Filename to save report information otherwise will print to stdout')

    parser_assess.add_argument('-s','--stats',action='store_true',
                               help='Only print the summary stats (L50,N50) and do not run completeness assay')

    parser_assess.add_argument('-m','--method',default="busco",choices=["busco"], # later FGMP
                               help='Method for assess completeness - BUSCO is default')
    
    parser_assess.add_argument('--clade',default='fungi_odb9',
                               help='BUSCO clade to use in completeness assessment')

    parser_assess.add_argument('--AAFTF_DB',type=str,
                               required=False,
                               help="Path to AAFTF resources, defaults to $AAFTF_DB")
    
    parser_assess.add_argument('--tmpdir',type=str,
                               required=False,default="working_AAFTF",
                               help="Temporary directory to store datafiles and processes in")

    parser.set_defaults(func=run_subtool)

    ### process args now ###
    # if no args then print help and exit
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()
    fh = logging.FileHandler('AAFTF.log')
    fh.setLevel(logging.INFO)
    
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    
    if args.quiet:
        logger.setLevel(logging.ERROR)
    try:
        args.func(parser, args)

    except IOError as e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
