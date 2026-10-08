# CLAUDE.md — Claude Code Project Context for AAFTF

See [AGENTS.md](AGENTS.md) for full development guidelines, code style, and common patterns.

## Key Architecture Reminders

- Entry point: `aaftf/main.py` (console script `AAFTF = "aaftf.main:main"`) — builds the top-level parser, wires up subcommands via `register_subcommands(parser)` (imported from `aaftf/_menu.py`; it creates the subparsers action and lists the subcommands in `AAFTF --help` under four group titles — Setup, Steps, Pipeline, Annotation — using the help-only `SubcommandGroup` action from `aaftf/utility.py`), sets up logging from `-q/--quiet` / `-v/--verbose` (`setup_logging()`), and dispatches by calling `args.func(**vars(args))`; `main()` returns a numeric exit code (0 ok, 1 error, 2 `FileNotFoundError`, 130 Ctrl-C) and subcommands raise exceptions instead of calling `sys.exit` (see AGENTS.md). There is no central dispatcher function: each subcommand parser binds its own subtool's `run` directly via `parser_x.set_defaults(func=<module>.run)` inside its `<name>_menu()` function in `_menu.py`. Subcommand parsers register no `aliases=[...]` — each subcommand has exactly one canonical name.
- All argparse subcommand-parser definitions ("menus") live in one place: `aaftf/_menu.py`. Each subcommand has a `<name>_menu(subparsers)` function there (e.g. `trim_menu`, `depth_menu`) that builds and registers that subtool's `subparsers.add_parser(...)` block and its arguments, then calls `parser_x.set_defaults(func=<module>.run)`; `aaftf/_menu.py` also defines `add_verbosity_args()` and `register_subcommands()` (which calls every `<name>_menu()` group by group, in `AAFTF --help` order, and returns the subparsers action), and imports every subcommand module to reference their `run` functions. Each subcommand module itself (`trim.py`, `depth.py`, etc.) only contains the `run(**kwargs)` execution logic — no parser-building code.
- Every subcommand's `run()` function takes its CLI dest names directly as keyword arguments (not an `args`/`parser` pair) and ends its parameter list with `**kwargs` to absorb argparse's bookkeeping attributes (`command`, `func`, and `quiet` where the step doesn't use it) that ride along in `vars(args)`. The body uses those parameters directly (e.g. `read1`, `workdir`); a parameter that needs to change during the run (e.g. `basename` auto-derived from `read1`, `workdir` defaulted when not given) is simply reassigned as a local variable. Only `assemble.py`'s internal helpers (`run_spades`, `run_megahit`, `run_unicycler`) take a forwarded subset of fields.
- `pipe` is a `run()`-only parameter (not a CLI option): the `pipeline_*` subcommands pass `pipe=True` to suppress each step's "next command" hint.
- Because `run()` never mutates a caller-supplied object, tests must not expect derived values (e.g. `trim.run()`'s auto-derived `basename`) to show up on the kwargs the test constructed — assert on side effects (subprocess commands, files written) instead. See `tests/test_trim.py::_run_bbduk` for the pattern.
- `aaftf/pipeline.py` has one entry point per pipeline subcommand (`run_short`, `run_long`, `run_hybrid`), each following the same convention and calling each step through `_run_step()` / `<module>.run(**_step_kwargs(name, shared, ...))`; shared stretches live in private helpers (`_trim_and_filter`, `_assemble_flye`, `_clean_and_finish`).
- `CustomHelpFormatter` lives in `aaftf/utility.py`. Every subcommand parser (built in `_menu.py`) uses `formatter_class=CustomHelpFormatter`.
- Working directories: `workdir, custom_workdir = make_workdir(workdir, "<name>")` at the start and `cleanup_workdir(workdir, debug, custom_workdir)` at the end (both in `utility.py`); `make_workdir` also opens the step's log file. The body after `make_workdir` is wrapped in `try:` … `finally: close_step_log()` (`utility.py`), so a step that raises still closes its log (the workdir is kept for inspection); `assemble` has no `make_workdir` and instead removes its auto-named assembler folder (`spades_*`, `megahit_*`, `unicycler_*`, `flye_*`) after the run unless `debug`. Logging is via a module-level `logger = logging.getLogger(__name__)`.

## CLI Framework Conventions

- Every subcommand parser has `-q/--quiet` (dest `quiet`) and `-v/--verbose` (dest `debug`, so code and `run()` parameters use `debug`). There is no `--pipe` CLI flag.
- `add_verbosity_args(target)` (in `aaftf/_menu.py`; formerly `menu_common_args`) adds exactly those two arguments to `target`. Every `<name>_menu()` function calls it **last**, passing its own `optional = parser_x.add_argument_group("optional arguments")` group (not the parser itself), so they render as the final entries of that subcommand's "optional arguments" section.
- All other args a subcommand needs (`-c/--cpus`, `-w/--workdir/--tmpdir`, `-1/--read1`/`-2/--read2`, etc.) are implemented locally within that subcommand's own `<name>_menu()` — there is no shared parent-parser mechanism (argparse `parents=[...]` shared `Action` objects across subcommands and per-subcommand overrides silently corrupted others). When adding such a block to a new menu function, copy it verbatim from a similar existing one (e.g. `vecscreen_menu`, `covpurge_menu`) to keep help text/defaults consistent. Path arguments use `type=str`.
- Each `<name>_menu()` function that has any required argument creates two argument groups — `required = parser_x.add_argument_group("required arguments")` and `optional = parser_x.add_argument_group("optional arguments")` — and adds each local argument to the appropriate one (required = `required=True` with no `default=`; everything else = optional).
- Adding a new subcommand: write its `<name>_menu(subparsers)` function in `aaftf/_menu.py` (create the subparser, add required/optional groups and any subtool-specific args, call `add_verbosity_args(optional)` last, then `parser_x.set_defaults(func=<module>.run)`), add it to the right group's list in `register_subcommands()`, and write the subtool's `run(**kwargs)` in its own module (named params for the dest names it uses, plus a trailing `**kwargs`).
- Testing a subcommand's CLI parsing without executing the tool: patch `aaftf.<module>.run` with a `side_effect` that captures `**kwargs` into a `Namespace`, then call `aaftf.main.main()` with `sys.argv` patched. See `tests/test_cli.py::_parse_with_main`.

## Subcommand Reference

"Module" below is where each subtool's `run(**kwargs)` execution logic lives; its parser is in `aaftf/_menu.py` (as `<name>_menu()`).

| Canonical name | Module | Key external tools |
|---|---|---|
| `dependency` | `aaftf/dependency.py` | (none — checks PATH for all tools below) |
| `database` | `aaftf/database.py` | (none — urllib only) |
| `trim` | `aaftf/trim.py` | bbduk.sh, trimmomatic, fastp |
| `mito` | `aaftf/mito.py` | NOVOPlasty.pl, minimap2, reformat.sh (subsampling) |
| `filter` | `aaftf/filter.py` | bbduk.sh, bowtie2, bwa, minimap2, samtools |
| `assemble` | `aaftf/assemble.py` | spades.py, megahit, unicycler, flye |
| `vecscreen` | `aaftf/vecscreen.py` | blastn, makeblastdb |
| `taxpurge` | `aaftf/taxpurge.py` | sourmash + sourmash_plugin_branchwater (`fastmultigather`, `tax genome`) |
| `covpurge` | `aaftf/covpurge.py` | bwa, samtools |
| `fcs_screen` | `aaftf/fcs_screen.py` | run_fcsadaptor.sh (singularity/docker) |
| `fcs_gx_purge` | `aaftf/fcs_gx_purge.py` | run_gx.py (NCBI FCS-GX) |
| `rmdup` | `aaftf/rmdup.py` | minimap2 |
| `polish` | `aaftf/polish.py` | polypolish+bwa; pypolca+bwa+samtools+freebayes; aaftf_polca.sh+bwa+samtools+freebayes; nextPolish2+minimap2+samtools+yak; racon+minimap2 |
| `sort` | `aaftf/sort.py` | (none — BioPython only) |
| `assess` | `aaftf/assess.py` | (none — BioPython only) |
| `depth` | `aaftf/depth.py` | samtools, mosdepth, minimap2 or bwa |
| `pipeline_short` | `aaftf/pipeline.py` (`run_short`) | tools of its steps |
| `pipeline_long` | `aaftf/pipeline.py` (`run_long`) | flye, racon + tools of its steps |
| `pipeline_hybrid` | `aaftf/pipeline.py` (`run_hybrid`) | flye/unicycler, racon, polypolish/pypolca + tools of its steps |
| `fix_tbl` | `aaftf/fix_tbl.py` | (none) |

## Module `run()` Parameters

Keyword parameters of each `run()` (besides the trailing `**kwargs`):

- **dependency**: `quiet`
- **database**: `databases`, `force`
- **trim**: `read1`, `read2`, `basename`, `method`, `cpus`, `memory`, `minlen`, `avgqual`, `trimmomatic_adaptors`, `trimmomatic_clip`, `trimmomatic_leadingwindow`, `trimmomatic_trailingwindow`, `trimmomatic_slidingwindow`, `trimmomatic_quality`, `merge`, `dedup`, `cutfront`, `cuttail`, `cutright`, `debug`, `pipe`
- **mito**: `read1`, `read2`, `out`, `workdir`, `minlen`, `maxlen`, `seed`, `subsample`, `memory`, `debug`, `pipe`
- **filter**: `read1`, `read2`, `workdir`, `cpus`, `screen_accessions`, `screen_urls`, `screen_local`, `basename`, `aligner`, `memory`, `debug`, `pipe`
- **assemble**: `out`, `method` (required, no default), `read1`, `workdir`, `cpus`, `memory`, `isolate`, `careful`, `assembler_args`, `tmpdir`, `read2`, `longreads`, `merged`, `longread_type`, `genome_size`, `debug`, `pipe`. `read1` is optional in the parser (required by spades/megahit/unicycler; `flye` needs `longreads` instead).
- **vecscreen**: `infile`, `outfile`, `workdir`, `cpus`, `percent_id`, `stringency`, `debug`, `pipe`
- **taxpurge**: `input`, `outfile`, `phylum`, `workdir`, `sourdb` (names in `resources.SOURMASH_DATABASES` — `gbk-fungi`/`gbk-fungi-lite`, `gtdb-reps`/`gtdb-reps-lite`, `human-hg38` — or paths of sketch collections / RocksDB indexes; None = `resources.TAXPURGE_DEFAULT_DATABASES`), `taxonomy_csv`, `kmer`, `min_hashes`, `cpus`, `taxonomy`, `debug`, `pipe`. Each database is searched separately at its own `scaled` (threshold = `min_hashes` x scaled bp; memory of the largest one) and each contig takes the classification with the most matched bp. `AAFTF database` downloads a sourmash database together with its lineage CSV (`with` in `DATABASES`); `AAFTF database taxpurge` gets the default set
- **covpurge**: `input`, `outfile`, `read1`, `read2`, `workdir`, `cpus`, `mincovpct`, `debug`, `pipe`
- **fcs_screen**: `infile`, `outfile`, `container_engine`, `workdir`, `image`, `prok`, `fcs_script`, `debug`, `pipe`
- **fcs_gx_purge**: `input`, `outfile`, `workdir`, `taxid`, `db`, `debug`, `pipe`
- **rmdup**: `input`, `out`, `workdir`, `cpus`, `percent_id`, `percent_cov`, `minlen`, `exhaustive`, `debug`, `pipe`
- **polish**: `infile`, `outfile`, `method`, `memory`, `cpus`, `read1`, `read2`, `longreads`, `workdir`, `debug`, `pipe`. Each `--method` (`polypolish` (default), `pypolca`, `polca` (bundled `aaftf_polca.sh`, a patched MaSuRCA `polca.sh` from `patches/polca/`, installed into `bin/` via hatch `shared-scripts`), `nextpolish2`, `racon`) has its own `run_<method>()` function in `aaftf/polish.py`.
- **sort**: `input`, `out`, `minlen`, `name`
- **assess**: `input`, `report`, `telomere_monomer`, `telomere_n_repeat`, `telomere_window`
- **depth**: `input`, `out`, `read1`, `read2`, `longreads`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `pipe`, `min_contig_len`, `no_plot`, `plot_format`
- **pipeline_short** (`run_short`): `read1`, `read2`, `basename`, `phylum`, `cpus`, `tmpdir`, `assembler_args`, `method`, `memory`, `minlen`, `screen_accessions`, `screen_urls`, `mincontiglen`, `workdir`, `sourdb`, `mincovpct`, `debug`, `quiet`
- **pipeline_long** (`run_long`): `longreads`, `basename`, `longread_type`, `genome_size`, `phylum`, `sourdb`, `cpus`, `assembler_args`, `memory`, `mincontiglen`, `workdir`, `debug`, `quiet`
- **pipeline_hybrid** (`run_hybrid`): `read1`, `longreads`, `basename`, `phylum`, `read2`, `longread_type`, `genome_size`, `method` (`flye`/`unicycler`), `cpus`, `tmpdir`, `assembler_args`, `memory`, `minlen`, `screen_accessions`, `screen_urls`, `mincontiglen`, `workdir`, `sourdb`, `mincovpct`, `debug`, `quiet`
- **fix_tbl**: `table`, `report`, `output`

## Pipeline step defaults

`pipeline_short` runs trim → filter → assemble → vecscreen → [taxpurge] → covpurge → rmdup → sort → assess (mito, fcs_screen, fcs_gx_purge, polish and depth are optional and not run).
`pipeline_long` runs assemble (flye) → vecscreen → [taxpurge] → rmdup → polish (racon) → sort → assess (no covpurge: it needs Illumina reads). In every pipeline `-p/--phylum` is optional and taxpurge runs only when it is given (`_taxpurge_options()`), so the pipelines work when the taxonomy is unknown.
`pipeline_hybrid` runs trim → filter → assemble (flye on the long reads) → vecscreen → [taxpurge] → covpurge → rmdup → polish (racon) → polish (polypolish for paired, pypolca for single-end filtered Illumina reads) → sort → assess; all pipelines clean (vecscreen/taxpurge/covpurge/rmdup) before polishing, via the shared `_clean_and_finish(assembly, ..., taxpurge_options, covpurge_options, polish_steps)`, so discarded contigs are never polished; with `--method unicycler` Unicycler gets both read types and both polish steps are skipped. A step whose output file already exists is skipped (`_run_step`). Each step is called with `_step_kwargs(name, shared, **step_options)`:
- it starts from that step's CLI defaults, read from its `_menu.py` parser (`_subcommand_defaults()`), so every `run()` parameter is present and the pipeline never drifts from `AAFTF <step>` defaults — do not hard-code step settings in `pipeline.py`
- `shared` pipeline options (`cpus`, `memory`, `workdir`, `debug`, `quiet`) go only to steps that have that option; a `None` value (e.g. no `--memory`) keeps the step default. Exception: `assemble` gets `{workdir}/assemble_{method}` (`_assemble_workdir()`), since SPAdes restarts from and MEGAHIT refuses an existing folder, and `filter` has already created the pipeline workdir
- `step_options` are the pipeline options specific to that step plus the file names chaining the steps; `pipe=True` is always set
- `tests/test_pipeline.py` checks, for all three pipelines, that each step receives exactly its CLI options, with defaults kept unless a pipeline option overrides them

## Pipeline Workflow

Recommended order for full assembly QC:

```
trim → [mito optional, PE only] → filter → assemble → vecscreen
  → [fcs_screen / fcs_gx_purge optional] → [taxpurge, needs a phylum] → covpurge → rmdup
  → [polish optional, mainly with long reads] → sort → assess → [depth optional]
```

Files `AAFTF pipeline_short -o {base}` writes (`pipeline_long` / `pipeline_hybrid` also write `{base}.flye.fasta`, `{base}.racon.fasta` and, hybrid only, `{base}.polypolish.fasta` / `{base}.pypolca.fasta`, or `{base}.unicycler.fasta` with `--method unicycler`):

| Step | Input | Output |
|---|---|---|
| trim | raw FASTQ | `{base}_1P.fastq.gz`, `{base}_2P.fastq.gz` |
| filter | trimmed FASTQ | `{base}_filtered_1.fastq.gz`, `{base}_filtered_2.fastq.gz` |
| assemble | filtered FASTQ | `{base}.{method}.fasta` |
| vecscreen | assembled FASTA | `{base}.vecscreen.fasta` |
| taxpurge (with `-p`) | vecscreen FASTA | `{base}.taxpurge.fasta` |
| covpurge | taxpurge (or vecscreen) FASTA | `{base}.covpurge.fasta` |
| rmdup | covpurge FASTA | `{base}.rmdup.fasta` |
| sort | rmdup FASTA | `{base}.final.fasta` |
| assess | sorted FASTA | printed stats (optional report file) |

Run standalone, `polish` defaults to `<input prefix>.polished.fasta` and `depth` to `coverage_stats.txt`.

## External Tool Versions (minimum known-good)

- samtools >= 1.13 (pixi pins >= 1.24; needed for `sort --write-index` and `flagstat -O tsv`). No samtools version branching: aligner → sorted BAM pipelines use `utility.align_to_sorted_bam()`, which pipes SAM straight into `samtools sort` and writes the `.bai` index in the same step
- minimap2 >= 2.17
- mosdepth >= 0.3
- sourmash >= 4.8 with sourmash_plugin_branchwater >= 0.9.15 (`taxpurge`: `fastmultigather` + `tax genome`)
- SPAdes >= 3.15 for `assemble`
- NCBI FCS-adaptor 0.5.5 (hard-coded in `resources.py`)

## Recent Additions (2026-03-02)

### `depth` subtool — coverage analysis

**New file:** `aaftf/depth.py`; registered in `main.py`.

**Workflow:**
1. Counts reads in each input FASTQ
2. Maps Illumina reads with `minimap2 -ax sr` (default) or `bwa mem` (`--aligner bwa`); maps long reads with `minimap2 -ax map-ont` (default) / `map-pb` / `map-hifi` (`--longread_type`)
3. Merges BAMs when both read types are provided (`samtools merge`)
4. Runs `samtools flagstat` per read type, then `mosdepth` on the combined BAM
5. Parses `mosdepth.summary.txt` for per-contig mean depths and `mosdepth.global.dist.txt` for coverage breadth
6. Flags contigs with mean depth > assembly_mean + 3×SD as possible contaminants/organelles (uses population SD; contigs are the full assembly population, not a sample)

**Key `args` attributes:** `input`, `out`, `read1`, `read2`, `longreads`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `min_contig_len`, `no_plot`, `plot_format`

## Audit Fixes (2026-05-02)

### Bugs fixed

| File | Line | Issue | Fix |
|---|---|---|---|
| `filter.py` | 225, 251, 273 | `tempfiles[3]` IndexError — only `tempfiles[0]` existed | Replaced with `unsorted_bam` path variable |
| `polish.py` | 82 | `polish_log = None` inside loop before `open()` | Removed erroneous reset; branches already set it |
| `polish.py` | 71 | `i` undefined if `--iterations 0` | Added `sys.exit(1)` guard when `iterations < 1` |
| `trim.py` | 135 | Infinite loop: `dirname("/") == "/"` | Added `if new_path == findpath: break` guard |
| `vecscreen.py` | 288 | Stale `start`/`end` in else-branch for multi-hit contigs | Moved `start, end = sorted(...)` before the if/else |
| `vecscreen.py` | 262 | `global contigs_to_remove` leaked state across calls | Changed to local variable in `run()` |
| `sourpurge.py` | 113 | `"nomatch" in cols` tested list membership — fragile | Changed to `cols[1].strip() == "nomatch"` |
| `fcs_screen.py` | 31 | Missing DB fell back to string literal `"AAFTF_DB"` | Changed to `status()` + `sys.exit(1)` |

### CLI consistency fixes (main.py)

- Added `--pipe` and `-v/--debug` to `fix_tbl` parser
- Added `--pipe` to `pipeline` parser
- Removed dead `runall` stub from `run_subtool()`
- Fixed assembler help text: "Default: unicycler" → "Default: spades"

### pipeline.py refactor

- `sort_args` and `assess_args` now use `create_namespace()` (consistent with all other steps)
- Assembly-specific attributes (`merged`, `isolate`, `careful`) now passed via `required_args` before `create_namespace()` call

### Bioinformatics enhancements

- `assess.py`: `find_telomere()` window is now configurable (`--telomere_window`, default 200 bp); added N-gap count and total N bases to stats output
- `rmdup.py`: now computes and reports both N50 and N75; filtering threshold uses N75 (as before)
- `depth.py`: added `--min_contig_len` (default 500 bp) to exclude short contigs from outlier depth analysis

## Running Tests

```bash
# Unit tests only (no external tools required)
pixi run -e dev pytest tests/ -m unit -v

# All tests (integration tests need the external tools, e.g. the pixi `complete`/`dev` environment)
pixi run -e dev pytest tests/ -v

# Lint/format (ruff, codespell, pydocstyle via pre-commit)
pixi run -e dev pre-commit run --all-files

# Compile-check all modules
python -m py_compile aaftf/*.py

# Docs (Sphinx + Graphviz diagrams; output in docs/_build/html)
pixi run -e dev docs         # build once
pixi run -e dev docs-serve   # live reload on http://127.0.0.1:8000
```
