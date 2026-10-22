# Changelog

## Unreleased (Development)

Changes since `v0.7.0-beta.4`. This is a large refactor: many command-line options changed, so
check the **Breaking changes** before updating scripts.

### Breaking changes

 - **`pipeline` renamed to `pipeline_short`** (same options and behaviour, no alias for the old name).

 - **`sourpurge` split into `taxpurge` and `covpurge`** (no alias for the old name). `taxpurge` is the
   sourmash taxonomy screen (`-p/--phylum`, `--sourdb`, `--min_hashes`, `-k/--kmer`,
   `--just-show-taxonomy`); `covpurge` is the bwa low-coverage screen (`-1/--read1` required,
   `-2/--read2`, `-mc/--mincovpct`, `-c/--cpus`). Run `taxpurge` first so contaminants don't skew
   `covpurge`'s coverage baseline. The pipelines write `{base}.taxpurge.fasta` and
   `{base}.covpurge.fasta` instead of `{base}.sourpurge.fasta`.
 - **`-p/--phylum` is optional in every pipeline**: without it `taxpurge` is skipped, so the pipelines
   also run when the organism's taxonomy is unknown. `pipeline_short` and `pipeline_hybrid` always run
   `covpurge`; `pipeline_long` does not (it needs Illumina reads).
 - **sourmash and its databases are now optional**: since the pipelines need them only with
   `-p/--phylum`, `AAFTF dependency` lists sourmash with the optional tools, `AAFTF database required`
   no longer downloads a sourmash database (fetch them with `AAFTF database taxpurge`), and sourmash moved from the
   pixi `default` environment (and `environment.yml`) to the `optional` feature (`complete`).

 - **`assemble --method` is required**: there is no default assembler any more (it was spades); pass
   `--method spades` to keep the old behaviour. The pipelines still choose one for you:
   `pipeline_short` uses spades and `pipeline_long` / `pipeline_hybrid` use flye.

 - **Package renamed** from `AAFTF` to `aaftf` for imports (`import aaftf`, `from aaftf.utility import ...`).
   The PyPI distribution and the `AAFTF` command keep their names.
 - **Read options renamed**: `-l/--left` and `-r/--right` are now `-1/--read1` and `-2/--read2` in
   `trim`, `mito`, `filter`, `assemble`, `sourpurge`, `polish`, `depth` and `pipeline`
   (`-r` remains `--report` in `assess` and `fix_tbl`).
 - **Subcommands renamed**: `check_dependencies` → `dependency`, `download` → `database`. Subcommand
   aliases (e.g. `asm`, `stats`, `purge`, `gx`) were removed; each subcommand has one name.
 - **Common flags**: the version flag is now `-V/--version`; `-v` is now `-v/--verbose` (was
   `-v/--debug`). The `--pipe` option was removed from every subcommand: "next command" hints are
   shown unless `-q/--quiet` is given.
 - **Removed options**: `assemble --method dipspades` and `--haplocontigs` (dipSPAdes is no longer
   packaged by SPAdes); `mito --reference` and `mito --starting` (circular mitochondrial genomes are
   always rotated to start at the cob gene); `fix_tbl` no longer reads stdin / writes stdout via `-`.
 - **`pipeline` no longer runs `polish` or `mito`**: it runs trim → filter → assemble → vecscreen →
   [taxpurge] → covpurge → rmdup → sort → depth → assess, and each step now uses its own command-line defaults
   (e.g. `assemble` with `--isolate`, `polish` default `polypolish`) unless a pipeline option overrides
   them. `pipeline -m/--memory` is now passed to every step that has `-m` (was 75% of RAM).
 - **Databases are downloaded only by `AAFTF database`**: `filter`, `vecscreen`, `taxpurge` and
   `fcs_screen` read them from the database folder (`filter` and `vecscreen` read `contamdb`, which
   `AAFTF database contamdb` builds) and stop with the `AAFTF database ...` command to
   run when one is missing (`filter -a/-u` still download what they are given). `$AAFTF_DB` may list
   several folders separated by `:`, and defaults to `~/.cache/aaftf`.
 - **`vecscreen` rewritten, following NCBI FCS-adaptor's rules**: one BLASTN search (NCBI's VecScreen
   settings) against the new `contamdb` database. Terminal hits (within 50 bp of an end) with a score
   >= 19, or >= 16 at the very end, are trimmed through to that end; internal hits are cut out only
   when Strong (score >= 30), so weak and moderate matches no longer split chromosomes (long-read
   contigs could be cut into hundreds of pieces); contigs more than 50% vector and pieces under
   200 bp are dropped. The Euk/Prok contamination and mitochondrial screens are removed, so
   `{prefix}.mitochondria.fasta` is no longer written, and so are `-pid/--percent_id` and
   `-s/--stringency`.
 - **`filter` and `vecscreen` share one contamination database, `contamdb`**: `contamdb.fa` (plus
   its BLAST database) holds UniVec_Core without its fungal/yeast-derived entries (names containing
   yeast, cerevisiae, pombe, pichia, kluyveromyces, aspergillus, neurospora, candida, 2-micron or
   shuttle), the PhiX genome and BBTools' `lambda.fa.gz` and `adapters.fa`, with duplicate and
   contained sequences removed. Build it with `AAFTF database contamdb` (part of `required`) before
   running either step. The `euks`, `proks` and `mitodb` databases are gone, `univec` is now
   UniVec_Core (was the full UniVec), and `filter` no longer screens BBTools'
   `sequencing_artifacts.fa.gz`.
 - **Pipeline output names**: `sort` now writes `{base}.sorted.fasta`; `{base}.final.fasta` is
   written by `depth`, with likely mitochondrial contigs tagged (see Added).
 - **Exit codes**: errors now exit with 1, a missing input file, tool or database with 2, and Ctrl-C
   with 130 (the traceback is shown with `-v`).
 - **Changed defaults / behaviour**:
   - `assess` default telomere monomer is now `TAAC{3,5}`.
   - `mito` subsamples to 1,500,000 read pairs by default (`--subsample 0` uses all reads).
   - `-m/--memory` is an integer number of GB everywhere (`assemble` used a string).
   - `vecscreen` writes its outputs next to `-o/--outfile` (was always the current directory).
   - `assess` reports BASES MASKED / PERCENT MASKED only when the assembly is soft-masked.
   - `fcs_gx_purge` drops only contigs FCS-GX marks EXCLUDE (TRIM/FIX/REVIEW rows are warned about
     and kept; previously every reported contig was dropped).

### Added

 - `pipeline_short`, `pipeline_long` and `pipeline_hybrid` can run NCBI FCS-adaptor: `--fcs_screen`
   (with `--container_engine`) runs `fcs_screen` instead of vecscreen (both screen for the same
   adaptors and vectors, so only one runs), writing `{base}.fcs_screen.fasta`. Off by default.
   NCBI FCS-GX (`fcs_gx_purge`) is not a pipeline step: its ~470 GB database must be in RAM, so run
   it separately on the finished assemblies.
 - **`assemble --method flye`**: Flye long-read assembly (`-lr/--longreads`, `--longread_type`
   default `nano-hq`, `--genome_size`); resumes a previous run in the workdir with `--resume`, and
   the next-command hint suggests `vecscreen` (polish after cleaning). `-1/--read1` is now required only for
   spades, megahit and unicycler. Flye is needed by `pipeline_long` / `pipeline_hybrid` (pixi `default`).
 - **`pipeline_long`**: long-read pipeline assemble (flye) → vecscreen → rmdup → polish (racon) →
   sort → depth → assess; with `-p/--phylum` taxpurge also runs before rmdup (no covpurge, which needs
   Illumina reads).
 - **`pipeline_hybrid`**: Illumina + long reads: trim → filter → assemble (flye) → vecscreen → [taxpurge] →
   covpurge → rmdup → polish (racon, then polypolish for paired or pypolca for single-end reads) → sort →
   depth → assess (contigs are cleaned before polishing, so discarded contigs are never polished); `--method unicycler` assembles both read types with Unicycler and skips polishing.

 - **`dependency`** (was `check_dependencies`): checks the tools and Python packages in two groups —
   what the default pipeline needs (an error if missing) and what only optional steps or non-default
   options use — printed in a fixed order with each tool's location or where a missing one is used;
   `-q` prints only the status. Suggests rebuilding bowtie2 for AVX2 when the conda build is installed.
 - **`database`** (was `download`): lists every reference database with its size, users and storage
   folder, split into databases the default pipeline needs and optional ones; downloads them by short
   name, file name, `required`, `optional` or `all`; `--force` re-downloads.
 - **`mito --subsample PAIRS`** (default 1,500,000; BBTools `reformat.sh`, fixed seed), and a
   "next command" hint that suggests `filter -s <mito.fasta>` to screen mitochondrial reads out.
 - **Per-step log files**: steps with a working directory write `<workdir>/<step>.log` (INFO and above
   even with `-q`, full tracebacks on errors), kept whenever the working directory is kept (own `-w`,
   `-v`, or a failed run); steps without one write `./<step>.log` with `-v`.
 - `-c/--cpus` and `-m/--memory` are lowered to what the machine/job can provide, with a warning.
 - Grouped `AAFTF -h` (Setup / Steps / Pipeline / Utilities) with optional steps marked "(Optional)".
 - `depth`: quantized coverage plots (heatmap and barplot, paginated by contig length), `--min_contig_len`,
   `--no-plot`, `--plot-format`; errors on samtools merge/index/flagstat failures.
 - pixi environments `default` (default pipeline), `complete` (every optional tool) and `dev` (AAFTF
   editable with the `dev` and `docs` extras), mirrored by `environment.yml` and `environment.dev.yml`;
   Python packages come from `[project]` dependencies (same layout as Phyling). A `docs` extra
   (`pip install ".[docs]"`) replaces `docs/requirements.txt`, and Read the Docs fetches the git tags
   so the docs show the real version; an `install-bowtie2` task that
   rebuilds bowtie2 with AVX2 support.
 - Documentation site (`docs/`, Read the Docs) with a page for every subcommand, including `dependency`.
 - GitHub Actions CI (`ci.yml`): pre-commit lint, unit tests on Python 3.10–3.14, integration tests in
   the locked pixi environment, and a strict Sphinx docs build.
 - conda package recipe (`ci/recipe/recipe.yaml`, rattler-build) and a publishing workflow.

 - **`AAFTF database contamdb`** builds the shared `filter`/`vecscreen` database (see Breaking
   changes). An existing `contamdb.fa` is kept, so you can add your own sequences to it, and
   `vecscreen` rebuilds its BLAST database whenever `contamdb.fa` is newer; `--force` rebuilds the
   FASTA from its sources. Once built, the downloaded UniVec_Core and PhiX are deleted (only from
   the folder AAFTF downloads into, never from another `$AAFTF_DB` folder) and listed as
   "fulfilled by contamdb".
 - Duplicate and contained sequences are removed with BBTools `dedupe.sh` (exact matches, both
   orientations, seed `k=10`) from `contamdb.fa`, and from `filter`'s combined database when
   `-a/-u/-s` add sequences; intermediate files are deleted.
 - **`depth` tags likely mitochondrial contigs** (nothing is removed): a contig is flagged MITO when
   its depth is >= `--mito_fold` (5) x the nuclear baseline (median depth of contigs >= N50), its GC
   is >= `--mito_gc_drop` (10) percentage points below the assembly's, and it is >= `--mito_min_len`
   (1,000) bp. The report gains a GC% column and a mitochondrial-candidates section,
   `<report>.mito_candidates.tsv` lists the evidence for every contig, and `--tagged FASTA` writes
   the assembly with `[location=mitochondrion]` (NCBI's source modifier) added to MITO contig
   headers.
 - Every pipeline runs `depth` after `sort` (Illumina reads in `pipeline_short`/`pipeline_hybrid`,
   long reads in `pipeline_long`) to write `{base}.final.fasta` with mitochondrial contigs tagged,
   plus `{base}.coverage_stats.txt` and `{base}.coverage_stats.mito_candidates.tsv`; if `depth`
   fails, a warning is logged and the sorted assembly is used untagged. `mosdepth` is therefore now one of the
   tools the default pipelines need (`AAFTF dependency`, `environment.yml`, pixi `default`).

### Changed

 - **`taxpurge` uses current sourmash databases and commands** instead of `sourmash lca
   classify` and its 2017 LCA databases. Contigs are searched with the branchwater plugin's
   `sourmash scripts fastmultigather`, one database at a time at its own resolution, and classified
   with `sourmash tax genome`; a contig matched in several databases takes the classification with
   the most matched sequence. `--sourdb` takes one or more of `gbk-fungi` / `gbk-fungi-lite`
   (GenBank fungi, `scaled=1000` / `10000`), `gtdb-reps` / `gtdb-reps-lite` (GTDB rs226
   representatives, bacteria and archaea) and `human-hg38`, or paths of sketch collections /
   RocksDB indexes (with `--taxonomy_csv`); it replaces `--sourdb_type`. The default,
   `gbk-fungi-lite gtdb-reps-lite human-hg38`, is ~0.7 GB (`AAFTF database taxpurge`) and needs
   ~1 GB of memory; the full-resolution databases (~9 GB) also classify shorter contigs.
   `--min_hashes` (default 2) sets the minimum match and `-c/--cpus` the threads. In a test (yeast
   assembly with *E. coli*, *C. acnes* and human contigs mixed in) the default removed 38 of 315
   bacterial and 33 of 40 human contigs and no yeast contig; the full databases removed 143 of 315
   bacterial contigs (the LCA version: 1 with its old defaults). The old `sm_gbk`, `sm_gtdbrep`
   and `sm_gtdb` LCA files are no longer used. Needs `sourmash_plugin_branchwater`.
 - `AAFTF database`: naming a sourmash database also downloads its lineage CSV, and `taxpurge`
   downloads `taxpurge`'s default databases.
 - Sourmash databases are installed as branchwater RocksDB indexes, which `taxpurge` searches from
   disk with much less memory (`gbk-fungi` from its index: ~1 GB): `AAFTF database gtdb-reps` /
   `gtdb-reps-lite` download and unpack the published index; `gbk-fungi`, `gbk-fungi-lite` and
   `human-hg38` are downloaded as sketch files, indexed locally and the sketch files removed
   (listed as "fulfilled by rocksdb"). `AAFTF database --no-build` keeps those sketch files
   instead, and `taxpurge` searches whichever is installed.
 - `AAFTF database` lists only the top level of the database folders under "Other files" (a
   subfolder as one entry), and lists each sketch file with its index on the row below it.
 - **Single-end reads in the pipelines**: `pipeline_short` and `pipeline_hybrid` stopped after
   `trim` on single-end reads, because they looked for `{base}_1P.fastq.gz` while `trim` writes
   `{base}_1U.fastq.gz`; they now use `_1U`. `filter` writes single-end output as
   `{base}_filtered_1.fastq.gz` with every aligner (was `_filtered_U` with bbduk and `_filtered`
   with the others), the same name as read 1 of a pair.
 - `filter` in a reused `-w` work directory rebuilds `contamdb.fa` when its sources change (e.g.
   other `-a/-u/-s`), and redoes the aligner index and alignment when they are older than it;
   before, an outdated database or BAM could be reused. The removed/kept read counts are now
   computed the same way for every aligner (from the filtered FASTQ).
 - `filter`: every `--aligner` now screens the same sequences, the shared `contamdb` database (which
   includes BBTools' lambda phage and adapter FASTAs), so bowtie2, bwa and minimap2 remove lambda
   too (before, only bbduk did, through its built-in `lambda` reference); bbduk no longer also
   loads BBTools' own PhiX or artifact references.
 - `polish` is optional and recommended only with long reads; `rmdup`'s next-command hint now suggests `sort`.
 - `fix_tbl` rewritten: features are clipped at trimmed contig ends and marked partial (`<`/`>`),
   features inside trimmed regions or on excluded contigs are dropped, minus-strand features and
   multi-interval features are handled, and the output is written via a temporary file so a failed
   run never leaves a truncated table.
 - `sourpurge`: the classification status must be exactly `found`; `--just-show-taxonomy` saves the
   classification CSV and cleans up its working directory.
 - `rmdup` reads the assembly once and matches contig IDs exactly.
 - `-h` output unified across subcommands: consistent help text and metavars (`FASTQ`, `FASTA`, `DIR`,
   `INT`, `GB`, `BP`, `PCT`, `PREFIX`) and a consistent option order (inputs, outputs, working
   directory, options, CPUs, memory, `-q`/`-v`, then tool-specific options).
 - Aligner → sorted BAM steps pipe straight into `samtools sort --write-index` (samtools >=1.13; no
   intermediate SAM files).
 - Bundled data (NOVOPlasty seed, config template, cob start sequence) is read with `importlib.resources`.
 - Containers build the pixi `complete` environment and report the correct version (`AAFTF_VERSION`
   build argument, converted to PEP 440).
 - Release workflow understands prerelease tags (e.g. `v0.7.0-beta.4` + `prerelease` → `v0.7.0-beta.5`);
   PyPI builds use the full git history so the published version matches the tag.

### Fixed

 - `AAFTF database` no longer downloads a file that is only needed to build a database already
   built: PhiX/UniVec_Core once `contamdb` exists, and a sourmash sketch file once its index exists
   (also with `--no-build`, which used to download it again); `--force` still downloads them.
 - `taxpurge --sourdb_type gtdb` / `gtdbrep` failed: they downloaded GTDB `.zip` collections,
   which `sourmash lca classify` cannot read. `taxpurge` now reads such collections directly.
 - `assemble --merged`: with paired reads SPAdes got the merged reads as a separate single-read
   library (`--s1`); they are now the paired library's merged reads (`--merged`), and single-end
   input uses `-s` (plus `--s 2` for merged reads) instead of the deprecated `--s1`/`--s2`/`--pe1-*`
   syntax. MEGAHIT ignored `--merged`; it is now passed with `-r` (comma-joined with single-end reads).
 - `assemble --method megahit` passed `-m` to MEGAHIT as bytes (32 GB became 32 bytes); it is now
   converted from GB. An existing MEGAHIT output folder is refused instead of failing inside MEGAHIT.
 - `assemble`: SPAdes restart passed `--restart-from last` as one argument; Unicycler dropped the paired
   reads when `--merged` was given; a missing assembly output now fails the step.
 - `trim`: the Trimmomatic adaptor search ignored `<prefix>/share/trimmomatic`, and the found adaptor
   file was never passed to Trimmomatic; missing reads or adaptors now raise errors.
 - `filter`: an unknown `--aligner` silently produced no output; the single-end hint suggested an
   `assemble` command that would not parse.
 - `fcs_screen` failed when FCS-adaptor had already created its output folder; FCS failures are reported.
 - `fcs_gx_purge` ignored `run_gx.py` failures.
 - `mito` picked NOVOPlasty's output by directory order (a circular assembly could be missed); a run
   with no assembly now fails.
 - `sourpurge` crashed when no N50 contig had coverage data (the filter is now skipped with a warning).
 - `assess` crashed on missing or empty input and never closed its report file.
 - `depth`: `-1` and `--longread_preset` were listed as required in `-h`; samtools failures were ignored.
 - `vecscreen` left input files open.
 - Trim/filter work around a BBDuk PairStreamer bug in paired mode; `pipeline` type mismatches that
   crashed `sourpurge` and `polish` were fixed.

### Removed

 - `assemble --method dipspades` and `--haplocontigs`; `mito --reference` and `--starting`; the `--pipe`
   option; subcommand aliases; the old `check_dependencies` / `download` names.
 - Unused utility functions and duplicated helpers across modules (consolidated in `aaftf/utility.py`).
 - Obsolete shell-script "tests", BUSCO/AUGUSTUS benchmark files and cluster-only test data from `tests/`.
 - Unused packages from the environments and conda recipe (pilon, masurca, diamond, flye, fastqc,
   hifiasm, kraken2, gfatools, taxonkit).

### Known issues

 - `taxpurge` classifies only about 1 in 7 contaminant contigs under 5 kb, and can misclassify a
   few short contigs of the target organism; short contaminant contigs are better caught by
   `covpurge` or `fcs_gx_purge`.

### Internal

 - Every `run()` takes its command-line option names as keyword arguments; all argparse menus live in
   `aaftf/_menu.py`; subcommands raise exceptions instead of calling `sys.exit`; logging uses the
   `logging` module.
 - Google-style docstrings and type hints throughout; mypy-clean; PEP 8 names; ruff targets Python 3.10.
 - Test suite rewritten and extended (about 700 tests, 96% line coverage), with tests isolated in
   temporary directories.

## 0.6.2

### Changed

 - Release with updated tools and checking for pipeline.

## 0.6.1 (tagged v0.6.1-alpha2)

### Tests

 - **New pytest test suite** (`tests/`): 144 unit tests, all passing.
   Requires `pip install pytest pytest-cov biopython`.  Run with `python -m pytest tests/`.
   - `tests/conftest.py` — shared fixtures (FASTA/FASTQ files, mosdepth summary data, argparse Namespaces)
   - `tests/test_utility.py` — pure-Python functions in `utility.py`: `myround`, `calcN50`, `RevComp`, `softwrap`, `checkfile`, `line_count`, `countfasta`, `fastastats`, `countfastq`, `which`, `SafeRemove`
   - `tests/test_assess.py` — `assess.py`: `revcomp`, `findTelomere` (fwd/rev/T2T/none), `genome_asm_stats` (GC%, N50, L50, telomere counts, gz input, file output), `run()`
   - `tests/test_sort.py` — `sort.run()`: descending-length order, header renaming, minlen filtering, sequence content preservation
   - `tests/test_fix_tbl.py` — `fix_tbl.py`: `parse_tbl` (blanks, comments, partial coords), `parse_adjustments` (left/right trim, multi-range), `fix_tbl` (trim_left shifts, trim_right clamps, pass-through)
   - `tests/test_depth.py` — `depth.py` pure functions: `count_fastq_reads`, `parse_mosdepth_summary`, `_coverage_breadth_from_dist`, outlier threshold arithmetic
   - `tests/test_cli.py` — CLI framework: `ALIAS_MAP` completeness, argparse defaults/required-args for `depth`/`assess`/`sort`, `--help` exits, alias routing
   - pytest configuration added to `pyproject.toml` (`[tool.pytest.ini_options]`); dev dependencies in `[project.optional-dependencies]`

### Added

 - **New `depth` subtool** (`AAFTF/depth.py`): calculates read depth of coverage for a genome assembly.
   Maps Illumina paired-end reads (via `minimap2 -ax sr` or `bwa mem`) and/or long reads (via `minimap2 -ax map-ont/map-pb/map-hifi`) to the assembly using a pipe chain (mapper stdout → `samtools sort`; no intermediate SAM files written to disk), then runs `mosdepth` to compute per-contig depth statistics.
   When both read types are provided their BAMs are merged before mosdepth.
   Writes `coverage_stats.txt` (configurable via `-o`) with three sections:
   1. Read input summary — per-file read counts and `samtools flagstat` alignment rates per read type
   2. Whole-assembly coverage — two mean depth estimates: mosdepth global (length-weighted) and per-contig arithmetic mean; plus % bases covered ≥1x
   3. Per-contig depth table sorted descending — two outlier tiers: `** OUTLIER` (mean > assembly mean + 3×SD, likely contaminant/organelle) and `ELEVATED` (2–3 SD above mean, worth inspecting)
   CLI: `AAFTF depth -i genome.sorted.fasta [-l fwd.fq] [-r rev.fq] [-lr longreads.fq] [-c cpus] [-o report]`
   Aliases: `coverage`, `cov`
   Required external tools: `minimap2` and/or `bwa`, `samtools`, `mosdepth`

### Fixed

 - **CLI alias routing**: Added `ALIAS_MAP` to `run_subtool()` in `AAFTF_main.py` so all argparse aliases (`asm`, `stats`, `dedup`, `pilon`, `polca`, `fix`, `purge`, `gx`, `trim_reads`, `read_trim`, `filter_reads`, `read_filter`, `vectorscreen`, `vector_blast`, `ncbi_fcs`, `ncbi_fcs-screen`, `ncbi_fcs-gx`, `ncbi_fcs_gx`, `mito_asm`, `mitochondria`) correctly dispatch to their canonical submodule instead of falling through to print-help-and-exit.
 - **`pipeline.py` assess step**: Fixed `AttributeError` — `assess_args` Namespace now includes `telomere_monomer` and `telomere_n_repeat` attributes required by `assess.run()`.
 - **`pipeline.py` assemble step**: Removed stale `asm_args.spades_tmpdir = None` assignment (wrong attribute name; `tmpdir` is already correctly populated via `assembleOpts`).
 - **Missing parser flags**: Added `-v/--debug` to `parser_assess`, `parser_sort`, and `parser_mito`; added `--pipe` to `parser_assess` and `parser_sort`; added `-v` shorthand to `parser_rmdup`'s `--debug` flag for consistency across all subcommands.
 - fcs_gx seemed like it was writing reports in the input fasta folder not the targeted output folder, should be fixed

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

## 0.4.1

### Fixed

 - Some bug fixes related to URLs.

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

## 0.2.4

### Added

 - conda/pypi packages.

## 0.2.3

### Added

 - Support a minimum length contig cutoff.

## 0.2.1

### Changed

 - Fix some README docs.
 - Sync zenodo with this release.

## 0.2.0

### Changed

 - Rework temporary file and prefix use throughout, not relying on a set working directory, so that multiple runs can be completed in a single folder without clashes of names.
 - Switch to BBTools bbduk for quality trimming instead of Trimmomatic.
 - Switch to BBTools bbduk instead of BWA/Bowtie for mapping to the contamination db - speed improvements and accuracy of k-mer matches.

### Added

 - Specify memory use for pilon and SPAdes assembly steps.
 - Support to generate scripts to run all steps in a pipeline with the 'pipeline' command.

## 0.1.0

 - Initial release, which includes automation for quality trimming and contaminant/adaptor removal, filtering, assembly, vector cleanup, refinement, duplication removal and sorting.
