# Changelog

## 0.6.1 (Development)

### Added

 - **New `depth` subtool** (`AAFTF/depth.py`): calculates read depth of coverage for a genome assembly.
   Maps Illumina paired-end reads (via `minimap2 -ax sr` or `bwa mem`) and/or long reads (via `minimap2 -ax map-ont/map-pb/map-hifi`) to the assembly, then runs `mosdepth` to compute per-contig depth statistics.
   When both read types are provided their BAMs are merged before mosdepth.
   Writes `coverage_stats.txt` (configurable via `-o`) with three sections:
   1. Read input summary — per-file read counts and `samtools flagstat` alignment rates per read type
   2. Whole-assembly coverage — mean depth and % bases covered (≥1x from mosdepth global distribution)
   3. Per-contig depth table sorted descending — contigs with mean depth > assembly_mean + 3×SD are flagged as `** OUTLIER (possible contaminant/organelle)`
   CLI: `AAFTF depth -i genome.sorted.fasta [-l fwd.fq] [-r rev.fq] [-lr longreads.fq] [-c cpus] [-o report]`
   Aliases: `coverage`, `cov`
   Required external tools: `minimap2` and/or `bwa`, `samtools`, `mosdepth`

### Fixed

 - **CLI alias routing**: Added `ALIAS_MAP` to `run_subtool()` in `AAFTF_main.py` so all argparse aliases (`asm`, `stats`, `dedup`, `pilon`, `polca`, `fix`, `purge`, `gx`, `trim_reads`, `read_trim`, `filter_reads`, `read_filter`, `vectorscreen`, `vector_blast`, `ncbi_fcs`, `ncbi_fcs-screen`, `ncbi_fcs-gx`, `ncbi_fcs_gx`, `mito_asm`, `mitochondria`) correctly dispatch to their canonical submodule instead of falling through to print-help-and-exit.
 - **`pipeline.py` assess step**: Fixed `AttributeError` — `assess_args` Namespace now includes `telomere_monomer` and `telomere_n_repeat` attributes required by `assess.run()`.
 - **`pipeline.py` assemble step**: Removed stale `asm_args.spades_tmpdir = None` assignment (wrong attribute name; `tmpdir` is already correctly populated via `assembleOpts`).
 - **Missing parser flags**: Added `-v/--debug` to `parser_assess`, `parser_sort`, and `parser_mito`; added `--pipe` to `parser_assess` and `parser_sort`; added `-v` shorthand to `parser_rmdup`'s `--debug` flag for consistency across all subcommands.

### Enhanced

 - **Version reporting now includes git checkout hash**: Enhanced version system to include short git commit hash (7 characters) for development installations. Version format examples:
   - Clean working tree: `0.6.0-alpha1-7-g261967e+261967e`
   - Dirty working tree: `0.6.0-alpha1-7-g261967e.dirty+261967e`
   - Tagged release: `0.6.0+261967e`
 - Version information is automatically detected from git repository for `pip install -e .` installations
 - Maintains backward compatibility with packaged installations and PEP 440 compliance
 - Improved development traceability and debugging support

## 0.6.0

## Added new features

 - unicycler as assembler for short reads
 - default spades mode is with --isolate (to turn off use --no-isolate --careful)
 - add menu item 'fix' / 'fix_tbl' which supports updating NCBI tbl format after fcs screening (have to run the `fcs clean genome` command)
## 0.5.0

### Added new features

 - NCBI fcs screening added for vector screening
 - NCBI fcs_gx screening for contamination (requires fast SSD disk or large memory)
 - Improve sourpurge to support database downloading, 3 types supported now (genbank 2017 microbial freeze, gtdb, gtdb_rep). Default for gtdb is rs214 release.
 - pilon tool renamed to polish and supports pilon, polca (masurca tool) is default, nextpolish. Racon not yet implemented for long read based polishing.
 - aliases for menu items (eg stats->assess; pilon->polish; asm->assemble)

### Bugs fixed

 - NCBI mitochondria genome download now points to the single FNA file instead of split between two

## 0.4.0

### Added

 - gzipped FastA files are supported by `AAFTF assess`
 - Telomere info reported in assess

## 0.3.3

### Added

 - pypi packaging for install

## 0.3.2

### Added

 - Support --careful and --isolate mode for spades run; default is --careful which seems to be a little better in tests

### Fixed

 - Add dependencies to requirements.txt and environment.yml

## 0.3.1

### Added
  - support for fastp in trimming (`trim`) command which can now support merging (--merge) and de-duplication (--dedup) as well as 5' and 3' trimming options of fastp (--cutfront, --cuttail, --cutright).
  -  added tests and summary stats to compare performance of different trimming/merging strategies on final genome and gene content

### Fixed
  - fix bug in `mito` to correctly spell novoplasty program that is used
  - revamped the resources.py to allow multiple files to represent mitochondria ref db for vector screening and alternative sourmash LCA databases - supporting gtdb-r207 and gtdb-rep-r207
  - spades in `assemble` command support merged long singleton reads along with paired end input


## 0.3.0

### Added

- Added NOVOplasty mediated mitochondrial assembly, `AAFTF mito`
- integrated `AAFTF mito` into `AAFTF pipeline`
- updated sourmash LCA database link and will download if not found


## 0.2.5

### Added

- Issue #10 added GC% in the assessment report table
- Added --mem option to pilon to up the heapsize for java runs
- merged changes in namespace by @gamcil and a tmpdir
