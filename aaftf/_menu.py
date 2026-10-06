"""Consolidated argparse subcommand parser (menu) definitions for AAFTF.

Each ``<name>_menu(subparsers)`` function below builds and registers one
AAFTF subcommand's argparse parser via ``subparsers.add_parser(...)``, and
binds that subtool's ``run()`` function to it via ``parser_x.set_defaults(
func=<module>.run)``. This keeps all CLI surface/wiring in one place,
separate from each subcommand module's own ``run(**kwargs)`` execution logic.
``aaftf.main`` invokes the selected subtool via ``args.func(**vars(args))``.

``add_verbosity_args()`` adds the two arguments common to every subcommand
(``-q/--quiet``, ``-v/--verbose``). Call it last, passing the
subcommand's own "optional arguments" group (not the parser itself), so these
common flags render as the final entries in that group instead of appearing
in a separate leading section.
"""

import argparse as ap

import aaftf.assemble as assemble
import aaftf.assess as assess
import aaftf.database as database
import aaftf.dependency as dependency
import aaftf.depth as depth
import aaftf.fcs_gx_purge as fcs_gx_purge
import aaftf.fcs_screen as fcs_screen
import aaftf.filter as aaftf_filter
import aaftf.fix_tbl as fix_tbl
import aaftf.mito as mito
import aaftf.pipeline as pipeline
import aaftf.polish as polish
import aaftf.rmdup as rmdup
import aaftf.sort as aaftf_sort
import aaftf.sourpurge as sourpurge
import aaftf.trim as trim
import aaftf.vecscreen as vecscreen
from aaftf.utility import CustomHelpFormatter, SubcommandGroup

__all__ = [
    "register_subcommands",
    "database_menu",
    "add_verbosity_args",
    "trim_menu",
    "mito_menu",
    "filter_menu",
    "assemble_menu",
    "vecscreen_menu",
    "fcs_screen_menu",
    "fcs_gx_purge_menu",
    "sourpurge_menu",
    "rmdup_menu",
    "polish_menu",
    "sort_menu",
    "assess_menu",
    "fix_tbl_menu",
    "depth_menu",
    "pipeline_short_menu",
    "pipeline_long_menu",
    "pipeline_hybrid_menu",
    "dependency_menu",
]


def register_subcommands(parser: ap.ArgumentParser) -> ap._SubParsersAction:
    """Add every AAFTF subcommand to ``parser``, listed in ``AAFTF --help`` under four group titles.

    Args:
        parser: The top-level ``AAFTF`` parser.

    Returns:
        The subparsers action holding all the subcommands.
    """
    subparsers = parser.add_subparsers(dest="command", metavar="<command>", help=ap.SUPPRESS, prog=parser.prog)
    groups = [
        ("Setup (dependencies and databases)", [dependency_menu, database_menu]),
        (
            "Steps",
            [trim_menu, mito_menu, filter_menu, assemble_menu, vecscreen_menu, sourpurge_menu, fcs_screen_menu, fcs_gx_purge_menu, rmdup_menu, polish_menu, sort_menu, assess_menu, depth_menu],
        ),
        ("Pipeline", [pipeline_short_menu, pipeline_long_menu, pipeline_hybrid_menu]),
        ("Annotation", [fix_tbl_menu]),
    ]
    for title, menus in groups:
        first = len(subparsers._choices_actions)
        for menu in menus:
            menu(subparsers)
        # argparse lists every subcommand in one block; show this group's entries under its own title instead
        parser.add_argument_group(title)._group_actions.append(SubcommandGroup(subparsers._choices_actions[first:]))
    return subparsers


def add_verbosity_args(target: ap._ActionsContainer) -> ap._ActionsContainer:
    """Add the -q/--quiet and -v/--verbose arguments that every AAFTF subcommand has.

    ``target`` is normally a subcommand's "optional arguments" group (from
    ``parser.add_argument_group("optional arguments")``). Call this only
    after all of that subcommand's own optional arguments have been added,
    so ``-q/--quiet`` and ``-v/--verbose`` render as the final
    entries of that group.

    Args:
        target: The argument group (or parser) to add the arguments to.

    Returns:
        ``target``, for chaining.
    """
    target.add_argument("-q", "--quiet", action="store_true", dest="quiet", help="Only show warnings and errors")
    target.add_argument("-v", "--verbose", action="store_true", dest="debug", help="Show debug messages and tool stderr, and keep temporary working directories")
    return target


def database_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the database subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``database`` subcommand parser.
    """
    parser_database = subparsers.add_parser(
        "database",
        formatter_class=CustomHelpFormatter,
        description=("List AAFTF reference databases (with no arguments), or download the named ones into the database folder ($AAFTF_DB, else ~/.cache/aaftf). Name databases by abbreviation or file name, or use 'all'."),
        help="List or download AAFTF reference databases",
    )
    parser_database.add_argument(
        "databases",
        nargs="*",
        metavar="DATABASE",
        help="Databases to download, by abbreviation (e.g. univec) or file name (e.g. UniVec), or 'all'. Omit to list the databases and where they are stored.",
    )
    optional = parser_database.add_argument_group("optional arguments")
    optional.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist",
    )
    add_verbosity_args(optional)
    parser_database.set_defaults(func=database.run)
    return parser_database


def trim_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the trim subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``trim`` subcommand parser.
    """
    parser_trim = subparsers.add_parser(
        "trim",
        formatter_class=CustomHelpFormatter,
        description="This command trims reads in FASTQ format to remove low quality reads and trim adaptor sequences",
        help="Trim FASTQ input reads",
    )

    required = parser_trim.add_argument_group("required arguments")
    optional = parser_trim.add_argument_group("optional arguments")

    required.add_argument(
        "-1",
        "--read1",
        type=str,
        metavar="FASTQ",
        required=True,
        help="Read 1 (forward) FASTQ, or single-end FASTQ",
    )

    optional.add_argument("-2", "--read2", type=str, metavar="FASTQ", help="Read 2 (reverse) FASTQ for paired-end data")

    optional.add_argument("-o", "--out", type=str, dest="basename", help="Output file prefix; default: read 1's file name up to its first '_' (or first '.' if it has none)", metavar="PREFIX")

    optional.add_argument("-ml", "--minlen", type=int, metavar="BP", default=75, help="Minimum read length to keep after trimming")

    optional.add_argument(
        "-aq",
        "--avgqual",
        type=int,
        metavar="INT",
        default=10,
        help="Average Quality of reads must be > than this",
    )

    optional.add_argument(
        "--cutfront",
        action="store_true",
        help="Run fastp 5' trimming based on quality. WARNING: this operation will interfere deduplication for SE data",
    )

    optional.add_argument(
        "--cuttail",
        action="store_true",
        help="Run fastp 3' trimming based on quality. WARNING: this operation will interfere deduplication for SE data",
    )

    optional.add_argument(
        "--cutright",
        action="store_true",
        help="Run fastp move a sliding window from front to tail, if meet one window with mean quality < threshold. \nWARNING: this operation will interfere deduplication for SE data",
    )

    optional.add_argument("--method", default="bbduk", choices=["bbduk", "trimmomatic", "fastp"], help="Trimming method", type=str)

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")
    optional.add_argument("-m", "--memory", type=int, dest="memory", default=8, help="Max memory in GB", metavar="GB")
    add_verbosity_args(optional)

    trimmomatic_group = parser_trim.add_argument_group(title="Trimmomatic options")

    trimmomatic_group.add_argument("--trimmomatic_adaptors", type=str, default="TruSeq3-PE.fa", help="Trimmomatic adaptor file")

    trimmomatic_group.add_argument("--trimmomatic_clip", type=str, default="2:30:10", help="Trimmomatic ILLUMINACLIP argument")

    trimmomatic_group.add_argument("--trimmomatic_leadingwindow", type=int, default=3, help="Trimmomatic window processing arguments")

    trimmomatic_group.add_argument("--trimmomatic_trailingwindow", type=int, default=3, help="Trimmomatic window processing arguments")

    trimmomatic_group.add_argument(
        "--trimmomatic_slidingwindow",
        type=str,
        default="4:15",
        help="Trimmomatic window processing arguments",
    )

    trimmomatic_group.add_argument("--trimmomatic_quality", default="phred33", help="Trimmomatic quality encoding -phred33 or phred64")

    fastp_group = parser_trim.add_argument_group(title="Fastp options")

    fastp_group.add_argument(
        "--dedup",
        action="store_true",
        help="Run fastp deuplication of fastq reads (default uses ~4gb mem)",
    )

    fastp_group.add_argument("--merge", action="store_true", help="Merge paired end reads")

    parser_trim.set_defaults(func=trim.run)
    return parser_trim


def mito_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the mito subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``mito`` subcommand parser.
    """
    parser_mito = subparsers.add_parser(
        "mito",
        description="De novo assembly of mitochondrial genome using NOVOplasty, takes PE Illumina adapter trimmed data.",
        help="(Optional) De novo assembly of mitochondrial genome",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_mito.add_argument_group("required arguments")
    optional = parser_mito.add_argument_group("optional arguments")

    required.add_argument("-1", "--read1", metavar="FASTQ", type=str, required=True, help="Read 1 (forward) FASTQ; mito needs paired-end reads")

    required.add_argument("-2", "--read2", metavar="FASTQ", type=str, required=True, help="Read 2 (reverse) FASTQ for paired-end data")

    optional.add_argument(
        "-s",
        "--seed",
        type=str,
        metavar="FASTA",
        help="Seed for NOVOPlasty: a mitochondrial sequence (e.g. a gene, or a related species' mitochondrial genome) the assembly is extended from. Default: the bundled Aspergillus nidulans cob fragment",
    )

    optional.add_argument("-o", "--out", type=str, default="mito.fasta", help="Output mitochondrial genome FASTA", metavar="FASTA")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--minlen", default=10000, type=int, help="Minimum expected genome size", metavar="BP")

    optional.add_argument("--maxlen", default=100000, type=int, help="Maximum expected genome size", metavar="BP")

    optional.add_argument(
        "--subsample",
        type=int,
        metavar="PAIRS",
        default=1_500_000,
        help="Randomly keep only this many read pairs (same pairs on every run) before assembling, since mitochondrial coverage is usually far higher than needed; 0 uses all reads",
    )

    optional.add_argument("-m", "--memory", type=int, dest="memory", default=8, help="Max memory in GB", metavar="GB")

    add_verbosity_args(optional)

    parser_mito.set_defaults(func=mito.run)
    return parser_mito


def filter_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the filter subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``filter`` subcommand parser.
    """
    parser_filter = subparsers.add_parser(
        "filter",
        description="Filter reads which match contaminant databases such as phiX",
        help="Filter contaminanting reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_filter.add_argument_group("required arguments")
    optional = parser_filter.add_argument_group("optional arguments")

    required.add_argument("-1", "--read1", metavar="FASTQ", type=str, required=True, help="Read 1 (forward) FASTQ, or single-end FASTQ")

    optional.add_argument("-2", "--read2", metavar="FASTQ", type=str, help="Read 2 (reverse) FASTQ for paired-end data")

    optional.add_argument("-o", "--out", dest="basename", type=str, help="Output file prefix; default: read 1's file name up to its first '_' (or first '.' if it has none)", metavar="PREFIX")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--aligner", default="bbduk", choices=["bbduk", "bowtie2", "bwa", "minimap2"], help="Aligner for mapping reads to the contaminant sequences", type=str)

    optional.add_argument("-a", "--screen_accessions", type=str, nargs="*", help="GenBank accession(s) whose sequences are screened out of the reads", metavar="ACCESSION")

    optional.add_argument("-u", "--screen_urls", type=str, nargs="*", help="URL(s) of FASTA files whose sequences are screened out of the reads", metavar="URL")

    optional.add_argument("-s", "--screen_local", type=str, nargs="+", help="Local FASTA file(s) to use contamination screen")

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    optional.add_argument("-m", "--memory", type=int, dest="memory", default=8, help="Max memory in GB", metavar="GB")

    add_verbosity_args(optional)

    parser_filter.set_defaults(func=aaftf_filter.run)
    return parser_filter


def assemble_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the assemble subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``assemble`` subcommand parser.
    """
    parser_asm = subparsers.add_parser(
        "assemble",
        description="Run assembler on cleaned reads: Illumina reads (spades, megahit), Illumina plus optional long reads (unicycler) or long reads only (flye)",
        help="Assemble reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_asm.add_argument_group("required arguments")

    required.add_argument("-o", "--out", type=str, required=True, help="Output assembly FASTA", metavar="FASTA")  # think about sensible replacement in future

    required.add_argument("--method", type=str, choices=["spades", "megahit", "unicycler", "flye"], required=True, help="Assembly method: spades or megahit (Illumina reads), unicycler (Illumina plus optional long reads) or flye (long reads)")

    shortread_group = parser_asm.add_argument_group(title="short-read (Illumina) arguments: spades, megahit, unicycler")

    shortread_group.add_argument("-1", "--read1", metavar="FASTQ", type=str, help="Read 1 (forward) FASTQ, or single-end FASTQ; required for spades, megahit and unicycler")

    shortread_group.add_argument("-2", "--read2", metavar="FASTQ", type=str, help="Read 2 (reverse) FASTQ for paired-end data")

    shortread_group.add_argument("--merged", type=str, dest="merged", help="Merged reads from flash or fastp or just single end reads")

    longread_group = parser_asm.add_argument_group(title="long-read (ONT/PacBio) arguments: flye, unicycler (hybrid)")

    longread_group.add_argument("-lr", "--longreads", type=str, help="Long-read FASTQ (PacBio or ONT); required for flye, optional for unicycler", metavar="FASTQ")

    longread_group.add_argument("--longread_type", choices=list(assemble.FLYE_READ_TYPES), default="nano-hq", help="Long-read type for flye")

    longread_group.add_argument("--genome_size", type=str, help="Estimated genome size for flye, e.g. 40m; default: estimated by flye", metavar="SIZE")

    optional = parser_asm.add_argument_group("optional arguments")

    optional.add_argument("-w", "--workdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--tmpdir", type=str, help="Temporary directory for the assembler", metavar="DIR")

    optional.add_argument("--assembler_args", action="append", help="Extra argument passed to the assembler (repeat for several)", metavar="ARG", type=str)

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    optional.add_argument("-m", "--memory", type=int, dest="memory", default=32, help="Max memory in GB for the assembler", metavar="GB")

    add_verbosity_args(optional)

    spades_group = parser_asm.add_argument_group(title="SPAdes options")

    spades_group.add_argument(
        "--no-careful",
        action="store_false",
        default=True,
        dest="careful",
        help="Disable --careful mode in spades (default: --careful is on)",
    )

    spades_group.add_argument(
        "--no-isolate",
        action="store_false",
        default=True,
        dest="isolate",
        help="Disable --isolate mode in spades (default: --isolate is on)",
    )

    parser_asm.set_defaults(func=assemble.run)
    return parser_asm


def vecscreen_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the vecscreen subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``vecscreen`` subcommand parser.
    """
    parser_vecscreen = subparsers.add_parser(
        "vecscreen",
        description="Screen contigs for vector and common contaminantion",
        help="BLASTN Vector and Contaminant Screening of contigs",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_vecscreen.add_argument_group("required arguments")
    optional = parser_vecscreen.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", type=str, required=True, dest="infile", help="Input genome assembly FASTA", metavar="FASTA")

    required.add_argument("-o", "--outfile", type=str, required=True, help="Output vector-screened assembly FASTA", metavar="FASTA")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("-pid", "--percent_id", type=int, help="Minimum percent identity for vector/contaminant BLAST hits", metavar="PCT")

    optional.add_argument("-s", "--stringency", default="high", choices=["high", "low"], help="Stringency to filter VecScreen hits")

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    add_verbosity_args(optional)

    parser_vecscreen.set_defaults(func=vecscreen.run)
    return parser_vecscreen


def fcs_screen_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the fcs_screen subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``fcs_screen`` subcommand parser.
    """
    parser_fcs_screen = subparsers.add_parser(
        "fcs_screen",
        description="Screen with NCBI fcs tool contigs for vector and common contaminantion",
        help="(Optional) NCBI Foreign Contaminant Screening for Vector sequences in contigs",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_fcs_screen.add_argument_group("required arguments")
    optional = parser_fcs_screen.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", type=str, required=True, dest="infile", help="Input genome assembly FASTA", metavar="FASTA")

    required.add_argument("-o", "--outfile", type=str, required=True, help="Output adaptor-screened assembly FASTA", metavar="FASTA")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--image", type=str, help="Container file (or will download and look in the database folder)")

    optional.add_argument(
        "--container_engine",
        type=str,
        default="singularity",
        choices=["singularity", "docker"],
        help="Container engine used to run fcs-adaptor",
    )

    optional.add_argument("--prok", action="store_true", help="Run in Prokaryote matching mode")

    optional.add_argument("--fcs_script", type=str, help="location of the run_fcsadaptor.sh script (or will download automatically)")

    add_verbosity_args(optional)

    parser_fcs_screen.set_defaults(func=fcs_screen.run)
    return parser_fcs_screen


def fcs_gx_purge_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the fcs_gx_purge subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``fcs_gx_purge`` subcommand parser.
    """
    parser_fcsgx = subparsers.add_parser(
        "fcs_gx_purge",
        description="Purge contigs based on fcs_gx results",
        help="(Optional) Purge contigs based on contamination search with fcs_gx",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_fcsgx.add_argument_group("required arguments")
    optional = parser_fcsgx.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", type=str, required=True, help="Input genome assembly FASTA", metavar="FASTA")

    required.add_argument("-o", "--outfile", type=str, required=True, help="Output FCS-GX-purged assembly FASTA", metavar="FASTA")  # think about sensible replacement in future

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument(
        "-t",
        "--taxid",
        type=int,
        default=4890,
        help="NCBI Taxonomy ID for contamaination matches, i.e. 4890 for Ascomycota",
    )

    optional.add_argument("-d", "--db", type=str, default="/my_tmpfs/gxdb/all", help="gxdb database path")

    add_verbosity_args(optional)

    parser_fcsgx.set_defaults(func=fcs_gx_purge.run)
    return parser_fcsgx


def sourpurge_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the sourpurge subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``sourpurge`` subcommand parser.
    """
    parser_sour = subparsers.add_parser(
        "sourpurge",
        description="Purge contigs based on sourmash results",
        help="Purge contigs based on sourmash results",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_sour.add_argument_group("required arguments")
    optional = parser_sour.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", type=str, required=True, help="Input genome assembly FASTA", metavar="FASTA")

    required.add_argument("-o", "--outfile", type=str, required=True, help="Output sourmash-purged assembly FASTA", metavar="FASTA")  # think about sensible replacement in future

    required.add_argument("-p", "--phylum", required=True, nargs="+", help="Phylum or phyla whose contigs are kept, e.g. Ascomycota", metavar="PHYLUM", type=str)

    optional.add_argument("-1", "--read1", metavar="FASTQ", type=str, help="Read 1 (forward) FASTQ, or single-end FASTQ")

    optional.add_argument("-2", "--read2", metavar="FASTQ", type=str, help="Read 2 (reverse) FASTQ for paired-end data")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--sourdb", type=str, help="sourmash LCA (k-31) taxonomy database; default: the one from 'AAFTF database'", metavar="FILE")

    optional.add_argument("-k", "--kmer", default="31", help="SourMash LCA kmersize when taxonomy database was built")

    optional.add_argument("-mc", "--mincovpct", default=5, type=int, help="Remove contigs whose coverage is below this percent of the N50 contigs' average coverage", metavar="PCT")

    optional.add_argument(
        "--sourdb_type",
        default="gbk",
        choices=["gbk", "gtdbrep", "gtdb"],
        help="Which sourpurge database to use.",
    )

    optional.add_argument("--just-show-taxonomy", dest="taxonomy", action="store_true", help="Show taxonomy information and exit")

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    add_verbosity_args(optional)

    parser_sour.set_defaults(func=sourpurge.run)
    return parser_sour


def rmdup_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the rmdup subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``rmdup`` subcommand parser.
    """
    parser_rmdup = subparsers.add_parser(
        "rmdup",
        description="Remove duplicate contigs",
        help="Remove duplicate contigs",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_rmdup.add_argument_group("required arguments")
    optional = parser_rmdup.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", type=str, required=True, help="Input genome assembly FASTA", metavar="FASTA")

    required.add_argument("-o", "--out", type=str, required=True, help="Output assembly FASTA with duplicate contigs removed", metavar="FASTA")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("-pid", "--percent_id", type=int, dest="percent_id", default=95, help="Minimum percent identity for a contig to count as a duplicate", metavar="PCT")

    optional.add_argument("-pcov", "--percent_cov", type=int, dest="percent_cov", default=95, help="Coverage of contig used to decide if it is redundant", metavar="PCT")

    optional.add_argument("-ml", "--minlen", type=int, default=500, help="Minimum contig length to keep", metavar="BP")

    optional.add_argument(
        "--exhaustive",
        action="store_true",
        help="Compute overlaps for every contig, otherwise only process contigs for L75 and below",
    )

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    add_verbosity_args(optional)

    parser_rmdup.set_defaults(func=rmdup.run)
    return parser_rmdup


def polish_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the polish subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``polish`` subcommand parser.
    """
    parser_polish = subparsers.add_parser(
        "polish",
        description="Polish contig sequences with Polypolish, pypolca, POLCA (polca.sh), NextPolish2 or Racon. Recommended for assemblies with long reads (or hybrid data); polishing a short-read-only assembly with the same short reads rarely helps and can introduce errors, so it is not part of AAFTF pipeline_short (pipeline_long and pipeline_hybrid run it).",
        help="(Optional) Polish contig sequences with short and/or long reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_polish.add_argument_group("required arguments")

    required.add_argument("-i", "--infile", "--input", type=str, dest="infile", required=True, help="Input genome assembly FASTA", metavar="FASTA")

    shortread_group = parser_polish.add_argument_group(title="polypolish / pypolca / polca / NextPolish2 required arguments")

    shortread_group.add_argument("-1", "--read1", metavar="FASTQ", type=str, help="Read 1 (forward) FASTQ, or single-end FASTQ; required for short-read polishing methods")

    shortread_group.add_argument("-2", "--read2", metavar="FASTQ", type=str, help="Read 2 (reverse) FASTQ for paired-end data; required for polypolish")

    longread_group = parser_polish.add_argument_group(title="NextPolish2 / Racon required arguments")

    longread_group.add_argument("-lr", "--longreads", type=str, help="Long-read FASTQ (PacBio or ONT); required for NextPolish2 (HiFi) and Racon", metavar="FASTQ")

    optional = parser_polish.add_argument_group("optional arguments")

    optional.add_argument("-o", "--out", "--outfile", type=str, dest="outfile", help="Output polished assembly FASTA", metavar="FASTA")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--method", type=str, choices=["polypolish", "pypolca", "polca", "nextpolish2", "racon"], default="polypolish", help="Polishing method")

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    add_verbosity_args(optional)

    pypolca_group = parser_polish.add_argument_group(title="pypolca / polca options")

    pypolca_group.add_argument("-m", "--memory", type=int, default=16, dest="memory", help="Max memory in GB", metavar="GB")

    parser_polish.set_defaults(func=polish.run)
    return parser_polish


def sort_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the sort subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``sort`` subcommand parser.
    """
    parser_sort = subparsers.add_parser(
        "sort",
        description="Sort contigs by length and rename FASTA headers",
        help="Sort contigs by length and rename FASTA headers",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_sort.add_argument_group("required arguments")
    optional = parser_sort.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", type=str, required=True, dest="input", help="Input genome assembly FASTA", metavar="FASTA")

    required.add_argument("-o", "--out", "--output", type=str, required=True, dest="out", help="Output sorted and renamed assembly FASTA", metavar="FASTA")

    optional.add_argument("-ml", "--minlen", type=int, default=0, help="Minimum contig length to keep", metavar="BP")

    optional.add_argument("-n", "--name", "--basename", type=str, default="scaffold", dest="name", help="Basename to rename FASTA headers")

    add_verbosity_args(optional)

    parser_sort.set_defaults(func=aaftf_sort.run)
    return parser_sort


def assess_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the assess subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``assess`` subcommand parser.
    """
    parser_assess = subparsers.add_parser(
        "assess",
        description="Assess completeness of genome assembly",
        help="Assess completeness of genome assembly",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_assess.add_argument_group("required arguments")
    optional = parser_assess.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", type=str, required=True, help="Input genome assembly FASTA", metavar="FASTA")

    optional.add_argument("-r", "--report", type=str, help="Write the report to this file (default: print to stdout)", metavar="FILE")

    optional.add_argument("-t", "--telomere_monomer", type=str, help="Telomere monomer pattern to search for.", default="TAAC{3,5}")

    optional.add_argument("-n", "--telomere_n_repeat", type=int, default=2, help="Telomere minimum number of monomer repeats.")

    optional.add_argument("--telomere_window", type=int, default=200, help="Number of bp to scan at each end for telomere repeats.")

    add_verbosity_args(optional)

    parser_assess.set_defaults(func=assess.run)
    return parser_assess


def fix_tbl_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the fix_tbl subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``fix_tbl`` subcommand parser.
    """
    parser_fix = subparsers.add_parser(
        "fix_tbl",
        description="Fix NCBI tbl file from a trim report from NCBI-FCS",
        help="Fix the TBL file offsets from trimmimg",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_fix.add_argument_group("required arguments")
    optional = parser_fix.add_argument_group("optional arguments")

    required.add_argument("-t", "--table", "--infile", type=str, required=True, help="Annotation table in NCBI .tbl format")

    required.add_argument("-r", "--report", type=str, required=True, help="NCBI FCS action report (tab-separated: accession, length, action, range(s), ...)", metavar="FILE")

    required.add_argument("-o", "--output", type=str, required=True, help="Write fixed TBL file")

    add_verbosity_args(optional)

    parser_fix.set_defaults(func=fix_tbl.run)
    return parser_fix


def depth_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the depth subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``depth`` subcommand parser.
    """
    parser_depth = subparsers.add_parser(
        "depth",
        description=("Calculate depth of coverage by mapping Illumina and/or long reads to a genome assembly with minimap2 (or bwa), then running mosdepth to compute per-contig depth statistics.  Contigs with mean depth > assembly_mean + 3*SD are flagged as possible contaminants or organellar sequences."),
        help="(Optional) Calculate read depth of coverage for genome assembly",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_depth.add_argument_group("required arguments")
    optional = parser_depth.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", type=str, required=True, dest="input", help="Input genome assembly FASTA", metavar="FASTA")

    optional.add_argument(
        "-1",
        "--read1",
        metavar="FASTQ",
        type=str,
        help="Read 1 (forward) FASTQ, or single-end FASTQ (Illumina); give -1 and/or -lr/--longreads",
    )

    optional.add_argument(
        "-2",
        "--read2",
        metavar="FASTQ",
        type=str,
        help="Read 2 (reverse) FASTQ for paired-end data (Illumina)",
    )

    optional.add_argument("-lr", "--longreads", type=str, help="Long-read FASTQ (PacBio or ONT)", metavar="FASTQ")

    optional.add_argument("-o", "--out", "--report", type=str, default="coverage_stats.txt", dest="out", help="Output coverage report", metavar="FILE")

    optional.add_argument("-w", "--workdir", "--tmpdir", type=str, dest="workdir", help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--aligner", default="minimap2", choices=["minimap2", "bwa"], help="Aligner for mapping Illumina reads", type=str)

    optional.add_argument(
        "--min_contig_len",
        type=int,
        default=500,
        metavar="BP",
        help="Minimum contig length to include in depth outlier analysis",
    )

    optional.add_argument(
        "--plot-format",
        choices=["pdf", "svg", "png"],
        default="pdf",
        dest="plot_format",
        help="Output format for coverage plots (requires matplotlib)",
    )

    optional.add_argument(
        "--no-plot",
        action="store_true",
        dest="no_plot",
        help="Disable coverage plot generation",
    )

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    add_verbosity_args(optional)

    longread_group = parser_depth.add_argument_group(title="long-read options")

    longread_group.add_argument(
        "--longread_preset",
        choices=["map-ont", "map-pb", "map-hifi"],
        dest="longread_preset",
        help="minimap2 preset for long reads; required with -lr/--longreads",
    )

    parser_depth.set_defaults(func=depth.run)
    return parser_depth


def pipeline_short_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the pipeline_short subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``pipeline_short`` subcommand parser.
    """
    parser_pipeline = subparsers.add_parser(
        "pipeline_short",
        description="Run the AAFTF short-read (Illumina) pipeline: trim, filter, assemble, vecscreen, sourpurge, rmdup, sort and assess. Each step uses its own defaults (see AAFTF <step> -h); only the options below override them.",
        help="Run the AAFTF pipeline on Illumina reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_pipeline.add_argument_group("required arguments")
    optional = parser_pipeline.add_argument_group("optional arguments")

    required.add_argument("-1", "--read1", metavar="FASTQ", type=str, required=True, help="Read 1 (forward) FASTQ, or single-end FASTQ")

    required.add_argument("-o", "--out", type=str, required=True, dest="basename", help="Output file prefix for every step's output files", metavar="PREFIX")

    required.add_argument("-p", "--phylum", required=True, nargs="+", help="Phylum or phyla whose contigs are kept, e.g. Ascomycota", metavar="PHYLUM", type=str)

    optional.add_argument("-2", "--read2", metavar="FASTQ", type=str, help="Read 2 (reverse) FASTQ for paired-end data")

    optional.add_argument("-w", "--workdir", type=str, help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--tmpdir", type=str, help="Temporary directory for the assembler", metavar="DIR")

    optional.add_argument("--assembler_args", action="append", help="Extra argument passed to the assembler (repeat for several)", metavar="ARG", type=str)

    optional.add_argument("--method", type=str, choices=["spades", "megahit", "unicycler"], default="spades", help="Assembly method")

    optional.add_argument("-ml", "--minlen", type=int, default=75, help="Minimum read length to keep after trimming", metavar="BP")

    optional.add_argument("-a", "--screen_accessions", type=str, nargs="*", help="GenBank accession(s) whose sequences are screened out of the reads", metavar="ACCESSION")

    optional.add_argument("-u", "--screen_urls", type=str, nargs="*", help="URL(s) of FASTA files whose sequences are screened out of the reads", metavar="URL")

    optional.add_argument("-mc", "--mincontiglen", type=int, default=500, help="Minimum length of contigs to keep", metavar="BP")

    optional.add_argument("--sourdb", type=str, help="sourmash LCA (k-31) taxonomy database; default: the one from 'AAFTF database'", metavar="FILE")

    optional.add_argument("--mincovpct", default=5, type=int, help="Remove contigs whose coverage is below this percent of the N50 contigs' average coverage", metavar="PCT")

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    optional.add_argument("-m", "--memory", type=int, dest="memory", help="Max memory in GB, passed to every step that has -m/--memory (trim, filter, assemble); default: each step's own default", metavar="GB")

    add_verbosity_args(optional)

    parser_pipeline.set_defaults(func=pipeline.run_short)
    return parser_pipeline


def pipeline_long_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the pipeline_long subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``pipeline_long`` subcommand parser.
    """
    parser_pipeline = subparsers.add_parser(
        "pipeline_long",
        description="Run the AAFTF long-read (ONT/PacBio) pipeline: assemble (flye), vecscreen, [sourpurge], rmdup, polish (racon), sort and assess. sourpurge runs only with -p/--phylum, as a taxonomy-only screen (its low-coverage filter needs Illumina reads). Each step uses its own defaults (see AAFTF <step> -h); only the options below override them.",
        help="Run the AAFTF pipeline on long reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_pipeline.add_argument_group("required arguments")
    optional = parser_pipeline.add_argument_group("optional arguments")

    required.add_argument("-lr", "--longreads", type=str, required=True, help="Long-read FASTQ (PacBio or ONT)", metavar="FASTQ")

    required.add_argument("-o", "--out", type=str, required=True, dest="basename", help="Output file prefix for every step's output files", metavar="PREFIX")

    optional.add_argument("--longread_type", choices=list(assemble.FLYE_READ_TYPES), default="nano-hq", help="Long-read type for flye")

    optional.add_argument("--genome_size", type=str, help="Estimated genome size for flye, e.g. 40m; default: estimated by flye", metavar="SIZE")

    optional.add_argument("-p", "--phylum", nargs="+", help="Phylum or phyla whose contigs are kept, e.g. Ascomycota; when given, sourpurge screens the contigs by taxonomy only (no coverage filter); default: sourpurge is not run", metavar="PHYLUM", type=str)

    optional.add_argument("--sourdb", type=str, help="sourmash LCA (k-31) taxonomy database; default: the one from 'AAFTF database'", metavar="FILE")

    optional.add_argument("-w", "--workdir", type=str, help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--assembler_args", action="append", help="Extra argument passed to the assembler (repeat for several)", metavar="ARG", type=str)

    optional.add_argument("-mc", "--mincontiglen", type=int, default=500, help="Minimum length of contigs to keep", metavar="BP")

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    optional.add_argument("-m", "--memory", type=int, dest="memory", help="Max memory in GB, passed to every step that has -m/--memory (assemble, polish); default: each step's own default", metavar="GB")

    add_verbosity_args(optional)

    parser_pipeline.set_defaults(func=pipeline.run_long)
    return parser_pipeline


def pipeline_hybrid_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the pipeline_hybrid subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``pipeline_hybrid`` subcommand parser.
    """
    parser_pipeline = subparsers.add_parser(
        "pipeline_hybrid",
        description="Run the AAFTF hybrid (Illumina + long-read) pipeline: trim, filter, assemble (flye on the long reads), vecscreen, sourpurge, rmdup, polish (racon, then polypolish for paired or pypolca for single-end Illumina reads), sort and assess. Polishing runs after the cleaning steps so discarded contigs are not polished. With --method unicycler both read types go to unicycler, which polishes its own assembly, so the polish steps are skipped. Each step uses its own defaults (see AAFTF <step> -h); only the options below override them.",
        help="Run the AAFTF pipeline on Illumina plus long reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_pipeline.add_argument_group("required arguments")
    optional = parser_pipeline.add_argument_group("optional arguments")

    required.add_argument("-1", "--read1", metavar="FASTQ", type=str, required=True, help="Read 1 (forward) FASTQ, or single-end FASTQ")

    required.add_argument("-lr", "--longreads", type=str, required=True, help="Long-read FASTQ (PacBio or ONT)", metavar="FASTQ")

    required.add_argument("-o", "--out", type=str, required=True, dest="basename", help="Output file prefix for every step's output files", metavar="PREFIX")

    required.add_argument("-p", "--phylum", required=True, nargs="+", help="Phylum or phyla whose contigs are kept, e.g. Ascomycota", metavar="PHYLUM", type=str)

    optional.add_argument("-2", "--read2", metavar="FASTQ", type=str, help="Read 2 (reverse) FASTQ for paired-end data")

    optional.add_argument("--method", type=str, choices=["flye", "unicycler"], default="flye", help="Assembly method")

    optional.add_argument("--longread_type", choices=list(assemble.FLYE_READ_TYPES), default="nano-hq", help="Long-read type for flye")

    optional.add_argument("--genome_size", type=str, help="Estimated genome size for flye, e.g. 40m; default: estimated by flye", metavar="SIZE")

    optional.add_argument("-w", "--workdir", type=str, help="Working directory for intermediate files; a temporary one is created and removed afterwards (kept with -v) when not given", metavar="DIR")

    optional.add_argument("--tmpdir", type=str, help="Temporary directory for the assembler", metavar="DIR")

    optional.add_argument("--assembler_args", action="append", help="Extra argument passed to the assembler (repeat for several)", metavar="ARG", type=str)

    optional.add_argument("-ml", "--minlen", type=int, default=75, help="Minimum read length to keep after trimming", metavar="BP")

    optional.add_argument("-a", "--screen_accessions", type=str, nargs="*", help="GenBank accession(s) whose sequences are screened out of the reads", metavar="ACCESSION")

    optional.add_argument("-u", "--screen_urls", type=str, nargs="*", help="URL(s) of FASTA files whose sequences are screened out of the reads", metavar="URL")

    optional.add_argument("-mc", "--mincontiglen", type=int, default=500, help="Minimum length of contigs to keep", metavar="BP")

    optional.add_argument("--sourdb", type=str, help="sourmash LCA (k-31) taxonomy database; default: the one from 'AAFTF database'", metavar="FILE")

    optional.add_argument("--mincovpct", default=5, type=int, help="Remove contigs whose coverage is below this percent of the N50 contigs' average coverage", metavar="PCT")

    optional.add_argument("-c", "--cpus", type=int, metavar="INT", default=1, help="Number of CPUs/threads to use")

    optional.add_argument("-m", "--memory", type=int, dest="memory", help="Max memory in GB, passed to every step that has -m/--memory (trim, filter, assemble, polish); default: each step's own default", metavar="GB")

    add_verbosity_args(optional)

    parser_pipeline.set_defaults(func=pipeline.run_hybrid)
    return parser_pipeline


def dependency_menu(subparsers: ap._SubParsersAction) -> ap.ArgumentParser:
    """Add the dependency subcommand parser.

    Args:
        subparsers: The top-level subparsers action to add the parser to.

    Returns:
        The new ``dependency`` subcommand parser.
    """
    parser_dependency = subparsers.add_parser(
        "dependency",
        formatter_class=CustomHelpFormatter,
        description="Check whether all external tool and Python package dependencies required by AAFTF are installed.",
        help="Check that AAFTF dependencies are installed",
    )

    optional = parser_dependency.add_argument_group("optional arguments")
    add_verbosity_args(optional)

    parser_dependency.set_defaults(func=dependency.run)
    return parser_dependency
