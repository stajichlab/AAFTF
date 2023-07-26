# Changelog

## 0.5.0

### Added new features

 - NCBI fcs screening added for vector screening
 - NCBI fcs_gx screening for contamination (requires fast SSD disk or large memory)
 - Improve sourpurge to support database downloading, 3 types supported now (genbank 2017 microbial freeze, gtdb, gtdb_rep). Default for gtdb is rs214 release.
 - pilon tool renamed to polish and supports pilon, polca (masurca tool), nextpolish. Racon not yet implemented for long read based polishing.
 - aliases for menu items (eg stats->assess; pilon->polisg; asm->assemble)

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
