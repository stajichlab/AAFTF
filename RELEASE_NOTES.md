Release notes for AAFTF 
=======================

# Automatic Assembly For The Fungi

* v0.3.1
  1. support for fastp in trimming (`trim`) command which can now support merging (--merge) and de-duplication (--dedup) as well as 5' and 3' trimming options of fastp (--cutfront, --cuttail, --cutright).
  2. fix bug in `mito` to correctly spell novoplasty program that is used
  3. revamped the resources.py to allow multiple files to represent mitochondria ref db for vector screening and alternative sourmash LCA databases - supporting gtdb-r207 and gtdb-rep-r207
  4. spades in `assemble` command support merged long singleton reads along with paired end input
  5. added tests and summary stats to compare performance of different trimming/merging strategies on final genome and gene content
  
* v0.3.0
  1. added novoplasty and `mito` sub-command
  2. better support for downloading sourmash lca db if it doesn't exist
  
* v0.2.4 
  1. conda/pypi packages

* v0.2.3
  1. support a minimum length contig cutoff

* v0.2.1
  1. Fix some README docs 
  2. Sync zenodo with this release

* v0.2.0
  1. rework temporary file, prefix use throughout, not relying on a set working directory so that multiple runs can be completed in a single folder
without clashes of names. 
  2. Switch to BBTools bbduk for quality trimming instead of Trimmomatic
  3. Switch to BBTools bbduk instead of BWA/Bowtie for mapping to contamination db. Speed improvements and accuracy of kmer matches to contamination over alignment-based cleanup.
  4. Specify memory use for pilon, spades assembly steps
  5. Support to generate scripts to run all steps in a pipeline with 'pipeline' cmd
  
* v0.1.0
  1. Initial release which includes automation for quality trimming and contaminant/adaptor removal, filtering, assembly, veccleanup, refinement, and duplication removal and sorting
