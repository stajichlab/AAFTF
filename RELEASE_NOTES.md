Release notes for AAFTF
=======================

# Automatic Assembly For The Fungi
* v0.6.2 (Stable)
  1. Release with updated tools and checking for pipeline.

* v0.6.1 (Development)
   **Tests**
   1. New pytest test suite (144 unit tests, all passing) in `tests/`
      - Requires: `pip install pytest pytest-cov biopython`
      - Run with: `python -m pytest tests/`
      - Includes unit tests for utility.py, assess.py, sort.py, fix_tbl.py, depth.py, and CLI framework

   **Added**
   1. New `depth` subtool - calculates read depth of coverage for genome assemblies
      - Maps Illumina paired-end reads via minimap2 (default) or bwa; long reads via minimap2 map-ont/map-pb/map-hifi
      - Mapper output is piped directly to `samtools sort` — no intermediate SAM files written to disk
      - Generates a coverage report with per-contig depth statistics and two outlier tiers:
        - `** OUTLIER`: mean depth > assembly mean + 3 SD (likely contaminant or organelle)
        - `ELEVATED`: mean depth 2–3 SD above mean (candidates worth inspecting)
      - Reports two mean depth estimates: mosdepth global (length-weighted) and per-contig arithmetic mean
      - Optional matplotlib coverage plots: per-scaffold heatmap, stacked bar chart, and depth histogram
      - Aliases: `coverage`, `cov`
      - Requires: minimap2 and/or bwa, samtools, mosdepth

   **Fixed**
   1. CLI alias routing - Fixed ALIAS_MAP in run_subtool() so all argparse aliases correctly dispatch to canonical submodules
   2. pipeline.py assess step - Fixed AttributeError by adding missing telomere_monomer and telomere_n_repeat attributes
   3. pipeline.py assemble step - Removed stale tmpdir assignment
   4. Missing parser flags - Added -v/--debug and --pipe flags for consistency across subcommands
   5. fcs_gx now writes reports to targeted output folder instead of input fasta folder

   **Enhanced**
   1. Version reporting now includes git checkout hash (7 characters) for development installations
      - Examples: 0.6.0-alpha1-7-g261967e+261967e (clean), 0.6.0-alpha1-7-g261967e.dirty+261967e (dirty)
      - Maintains PEP 440 compliance and backward compatibility

* v0.6.0
   1. add unicycler support
   2. default spades mode is with --isolate (use --no-isolate --careful to turn off)
   3. add 'fix' / 'fix_tbl' menu item for updating NCBI tbl format after fcs screening

* v0.5.0
   1. fcs NCBI tool supported with docker/singularity run
   2. fcs_gx from NCBI supported with docker/singularity - more accurate than sourpurge in many cases
   3. renamed pilon to polish step, defaults to POLCA (but pilon can still be used)
   4. Improve sourpurge to support database downloading - 3 types now supported (genbank 2017 microbial freeze, gtdb, gtdb_rep with rs214 default)
   5. Bug fix: NCBI mitochondria genome download now points to single FNA file instead of split between two

* v0.4.1
   1. Some bug fixes related to URLs

* v0.4.0
   1. linted code better with flake8
   2. Support gzipped assembly files in assess
   3. Add telomere counting to assess

* v0.3.3
   1. Created pypi packaging now

* v0.3.2
   1. In assemble; Support --isolate and --careful option; --careful is default (as was in previous versions)

* v0.3.1
   1. Support novoplasty Mitochondria targeted assembly runs with the `AAFTF mito` option (uses default A.nidulans seed or user-specified target)
   2. Support for alternative sourpurge databases - gtdb (r207) and gtdb-rep (r207) along with default 2018 genbank-k31 release
   3. trimming supports fastp now with support for de-duplication option; merged option to merge paired end overlapping reads; and additional novoseq/NextSeq trailing 'G' option trimming (--cutfront, --cuttail, --cutright)
   4. fix bug in `mito` to correctly spell novoplasty program that is used
   5. revamped the resources.py to allow multiple files to represent mitochondria ref db and alternative sourmash LCA databases - supporting gtdb-r207 and gtdb-rep-r207
   6. spades in `assemble` command support merged long singleton reads along with paired end input
   7. added tests and summary stats to compare performance of different trimming/merging strategies on final genome and gene content

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
   1. rework temporary file, prefix use throughout, not relying on a set working directory so that multiple runs can be completed in a single folder without clashes of names
   2. Switch to BBTools bbduk for quality trimming instead of Trimmomatic
   3. Switch to BBTools bbduk instead of BWA/Bowtie for mapping to contamination db - speed improvements and accuracy of kmer matches
   4. Specify memory use for pilon, spades assembly steps
   5. Support to generate scripts to run all steps in a pipeline with 'pipeline' cmd

* v0.1.0
   1. Initial release which includes automation for quality trimming and contaminant/adaptor removal, filtering, assembly, veccleanup, refinement, and duplication removal and sorting
