# AAFTF Development Guidelines

This document provides guidelines for coding agents (Claude Code, Codex, Copilot, ...) working on the AAFTF (Automatic Assembly For The Fungi) codebase. It is the single source of these instructions: `CLAUDE.md` only imports it.

## Project Overview

AAFTF is a Python-based bioinformatics toolkit for automated genome assembly, cleanup, mitochondrial genome assembly, and polishing focused on fungal genomes. The package uses external bioinformatics tools and provides a unified interface for genome assembly workflows.

## Build and Development Commands

### Installation and Setup
```bash
# pixi (recommended; locked versions in pixi.lock). Environments: default (pipeline tools),
# complete (+ optional tools; feature "optional"), dev (+ pytest, pre-commit, ruff; feature "dev")
pixi install -e dev
pixi run -e dev AAFTF -h
pixi run -e complete install-bowtie2   # optional: build bowtie2 from source (AVX2 issue)

# or conda: environment.yml (pipeline tools, AAFTF from PyPI) /
# environment.dev.yml (all tools, AAFTF installed editable from this checkout)
conda env create -f environment.dev.yml && conda activate aaftf-dev

# the package is built with hatchling + hatch-vcs; Python >=3.11,<3.14 (3.10: no sourmash_plugin_branchwater build; 3.14: no flye build)
pip install -e .
```

### Code Quality and Linting
```bash
# ruff settings are in ruff.toml (target py311)
# Run ruff for code style, linting, and import sorting (replaces flake8, isort, pyupgrade)
ruff check aaftf/
ruff format aaftf/

# Run ruff with automatic fixes
ruff check --fix aaftf/

# Run pydocstyle for documentation standards
pydocstyle --convention=google aaftf/

# Run codespell for spelling checks
codespell aaftf/ --ignore-words-list=nd,reacher,thist,ths,ure,referenc,wile,nin,pid

# All pre-commit hooks (recommended)
pre-commit run --all-files
```

### Testing
```bash
# Unit tests (no external tools needed)
pixi run -e dev pytest tests/ -m unit -v

# All tests (integration tests need the external tools, e.g. from the pixi environment)
pixi run -e dev pytest tests/ -v
```

### Docs
```bash
# Sphinx + Graphviz diagrams; output in docs/_build/html
pixi run -e dev docs         # build once
pixi run -e dev docs-serve   # live reload on http://127.0.0.1:8000
```

### Compile check
```bash
python -m py_compile aaftf/*.py
```

## Architecture

- Entry point: `aaftf/main.py` (console script `AAFTF = "aaftf.main:main"`) — builds the top-level parser, wires up subcommands via `register_subcommands(parser)` (imported from `aaftf/_menu.py`; it creates the subparsers action and lists the subcommands in `AAFTF --help` under four group titles — Setup, Steps, Pipeline, Utilities — using the help-only `SubcommandGroup` action from `aaftf/utility.py`), sets up logging from `-q/--quiet` / `-v/--verbose` (`setup_logging()`), and dispatches by calling `args.func(**vars(args))`; `main()` returns a numeric exit code (0 ok, 1 error, 2 `FileNotFoundError`, 130 Ctrl-C) and subcommands raise exceptions instead of calling `sys.exit` (see "Error Handling" below). There is no central dispatcher function: each subcommand parser binds its own subtool's `run` directly via `parser_x.set_defaults(func=<module>.run)` inside its `<name>_menu()` function in `_menu.py`. Subcommand parsers register no `aliases=[...]` — each subcommand has exactly one canonical name.
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
- **database**: `databases`, `force`, `no_build`
- **trim**: `read1`, `read2`, `basename`, `method`, `cpus`, `memory`, `minlen`, `avgqual`, `trimmomatic_adaptors`, `trimmomatic_clip`, `trimmomatic_leadingwindow`, `trimmomatic_trailingwindow`, `trimmomatic_slidingwindow`, `trimmomatic_quality`, `merge`, `dedup`, `cutfront`, `cuttail`, `cutright`, `debug`, `pipe`
- **mito**: `read1`, `read2`, `out`, `workdir`, `minlen`, `maxlen`, `seed`, `subsample`, `memory`, `debug`, `pipe`
- **filter**: `read1`, `read2`, `workdir`, `cpus`, `screen_accessions`, `screen_urls`, `screen_local`, `basename`, `aligner`, `memory`, `debug`, `pipe`
- **assemble**: `out`, `method` (required, no default), `read1`, `workdir`, `cpus`, `memory`, `isolate`, `careful`, `assembler_args`, `tmpdir`, `read2`, `longreads`, `merged`, `longread_type`, `genome_size`, `debug`, `pipe`. `read1` is optional in the parser (required by spades/megahit/unicycler; `flye` needs `longreads` instead).
- **vecscreen**: `infile`, `outfile`, `workdir`, `cpus`, `debug`, `pipe`. One blastn search with the VecScreen preset (`VECSCREEN_TASK_ARGS`; blastn has no `-task vecscreen`) against the `contamdb` BLAST database (`database.require_contamdb(blast=True)`, which re-runs `makeblastdb` when the index is missing or older than `contamdb.fa`), acted on with FCS-adaptor's rules: terminal hits (within 50 bp of an end) with score >= 19 (>= 16 at the tip) trimmed to the end, internal hits only when Strong (>= 30), contigs > 50% vector and pieces < 200 bp dropped. No Euk/Prok/mito screens. `AAFTF database contamdb` (`database.build_contamdb()`) builds `contamdb.fa` = NCBI `adaptors_for_screening_euks.fa` (`euk_adaptors`) + PhiX + BBTools `lambda.fa.gz`/`adapters.fa` (`utility.bbtools_resources()`), deduplicated with `dedupe.sh`, then `makeblastdb`; no UniVec/cloning vectors, whose fungal-derived parts (ADH1, LEU2, 2-micron, CEN/ARS) split fungal contigs; an existing `contamdb.fa` is kept (users may add sequences), `--force` rebuilds it. `filter` uses the same `contamdb.fa` (`require_contamdb()`; only that file is checked) as the base of its workdir `contamdb.fa`, plus `-a`/`-u`/`-s` extras
- **taxpurge**: `input`, `outfile`, `phylum`, `workdir`, `sourdb` (names in `resources.SOURMASH_DATABASES` — `gbk-fungi`/`gbk-fungi-lite`, `gtdb-reps`/`gtdb-reps-lite`, `human-hg38` — or paths of sketch collections / RocksDB indexes; None = `resources.TAXPURGE_DEFAULT_DATABASES`), `taxonomy_csv`, `kmer`, `min_hashes`, `cpus`, `taxonomy`, `debug`, `pipe`. Each database is searched separately at its own `scaled` (threshold = `min_hashes` x scaled bp; memory of the largest one) and each contig takes the classification with the most matched bp. A named database's RocksDB index (`rocksdb` in `SOURMASH_DATABASES`, found with `utility.find_db_index()`, which requires its `CURRENT` file) is used instead of its sketch file when present. `AAFTF database NAME` installs the index: `database._install_index_archive()` downloads and unpacks the published one (`rocksdb_url`; the DATABASES entry is then the index folder, `index`), or `database._build_index()` builds it from the sketch file (`build_index`) and removes the sketch file; `--no-build` (`no_build`) keeps the sketch file instead. A sourmash database comes with its lineage CSV (`with` in `DATABASES`); `AAFTF database taxpurge` gets the default set
- **covpurge**: `input`, `outfile`, `read1`, `read2`, `workdir`, `cpus`, `mincovpct`, `debug`, `pipe`
- **fcs_screen**: `infile`, `outfile`, `container_engine`, `workdir`, `image`, `prok`, `fcs_script`, `debug`, `pipe`
- **fcs_gx_purge**: `input` (list of assemblies, `-i` takes several), `db` (required; prefix, `.gxi` or folder, `_resolve_db()`), `taxid` (required), `outdir`, `workdir`, `cpus` (None = `available_cpus()`, passed as `GX_NUM_CORES`), `copy_to` (None = `$TMPDIR`, else `FALLBACK_COPY_TO` = `/scratch`), `debug`. A standalone utility ("Utilities" group), not a pipeline step: before anything runs it checks `run_gx.py` is on PATH and `get_ram()` (respects SLURM/cgroup limits) against the database size + `MEMORY_HEADROOM` (raises `MemoryError`; skipped when the db is on tmpfs, `_in_ram()`). A db on tmpfs, or one whose largest file passes `read_latency_ms()` (mean random 4 KiB read <= `FAST_READ_MS`), is used in place; otherwise all `{db}.*` files are copied once to a per-job `aaftf-gxdb-<SLURM_JOB_ID>-*` folder in `copy_to` (skipped when `copy_to` is on the db's filesystem, warned when on a `NETWORK_FS`; `_copy_target_ok()`), used for every assembly and removed in `finally` (no SIGTERM handling or stale-copy cleanup by design: use a `copy_to` the cluster removes per job). Writes `<prefix>.fcsgx_cleaned.fasta` (`cleaned_name()`) and `<prefix>.fcs_gx_report.txt` per input, next to it or in `outdir`
- **rmdup**: `input`, `out`, `workdir`, `cpus`, `percent_id`, `percent_cov`, `minlen`, `exhaustive`, `debug`, `pipe`
- **polish**: `infile`, `outfile`, `method`, `memory`, `cpus`, `read1`, `read2`, `longreads`, `workdir`, `debug`, `pipe`. Each `--method` (`polypolish` (default), `pypolca`, `polca` (bundled `aaftf_polca.sh`, a patched MaSuRCA `polca.sh` from `patches/polca/`, installed into `bin/` via hatch `shared-scripts`), `nextpolish2`, `racon`) has its own `run_<method>()` function in `aaftf/polish.py`.
- **sort**: `input`, `out`, `minlen`, `name`
- **assess**: `input`, `report`, `telomere_monomer`, `telomere_n_repeat`, `telomere_window`
- **depth**: `input`, `out`, `read1`, `read2`, `longreads`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `pipe`, `min_contig_len`, `no_plot`, `plot_format`, `tagged`, `mito_fold`, `mito_gc_drop`, `mito_min_len`. Flags MITO contigs (depth >= `mito_fold` x median depth of contigs >= N50, GC <= assembly GC - `mito_gc_drop`, length >= `mito_min_len`; `_mito_candidates()`), writes `<report>.mito_candidates.tsv`, and with `tagged` writes the assembly with `[location=mitochondrion]` (`MITO_TAG`) appended to MITO contig headers; nothing is removed
- **pipeline_short** (`run_short`): `read1`, `read2`, `basename`, `phylum`, `cpus`, `tmpdir`, `assembler_args`, `method`, `memory`, `minlen`, `screen_accessions`, `screen_urls`, `mincontiglen`, `workdir`, `sourdb`, `fcs_screen_run`, `container_engine`, `mincovpct`, `debug`, `quiet`
- **pipeline_long** (`run_long`): `longreads`, `basename`, `longread_type`, `genome_size`, `phylum`, `sourdb`, `fcs_screen_run`, `container_engine`, `cpus`, `assembler_args`, `memory`, `mincontiglen`, `workdir`, `debug`, `quiet`
- **pipeline_hybrid** (`run_hybrid`): `read1`, `longreads`, `basename`, `phylum`, `read2`, `longread_type`, `genome_size`, `method` (`flye`/`unicycler`), `cpus`, `tmpdir`, `assembler_args`, `memory`, `minlen`, `screen_accessions`, `screen_urls`, `mincontiglen`, `workdir`, `sourdb`, `fcs_screen_run`, `container_engine`, `mincovpct`, `debug`, `quiet`
- **fix_tbl**: `table`, `report`, `output`

## Pipeline step defaults

`pipeline_short` runs trim → filter → assemble → vecscreen → [taxpurge] → covpurge → rmdup → sort → depth → assess (mito, fcs_screen and polish are optional and not run; fcs_gx_purge is never run by a pipeline).
`pipeline_long` runs assemble (flye) → vecscreen → [taxpurge] → rmdup → polish (racon) → sort → depth (long reads) → assess (no covpurge: it needs Illumina reads). In every pipeline `-p/--phylum` is optional and taxpurge runs only when it is given (`_taxpurge_options()`), so the pipelines work when the taxonomy is unknown. Every pipeline can run `fcs_screen` (`--fcs_screen`, dest `fcs_screen_run`) instead of vecscreen (`_fcs_options()`; both screen for the same adaptors/vectors, so only one runs), writing `{base}.fcs_screen.fasta` in place of `{base}.vecscreen.fasta`. `fcs_gx_purge` is not a pipeline step: its ~470 GB database must be loaded into RAM, so it is run separately on finished assemblies.
`pipeline_hybrid` runs trim → filter → assemble (flye on the long reads) → vecscreen → [taxpurge] → covpurge → rmdup → polish (racon) → polish (polypolish for paired, pypolca for single-end filtered Illumina reads) → sort → depth (Illumina reads) → assess; `depth` writes `{base}.final.fasta` from `{base}.sorted.fasta` with likely mitochondrial contigs tagged (`_tag_mito()`; on failure it only warns and copies the sorted file); all pipelines clean (vecscreen/taxpurge/covpurge/rmdup) before polishing, via the shared `_clean_and_finish(assembly, ..., taxpurge_options, covpurge_options, polish_steps)`, so discarded contigs are never polished; with `--method unicycler` Unicycler gets both read types and both polish steps are skipped. A step whose output file already exists is skipped (`_run_step`). Each step is called with `_step_kwargs(name, shared, **step_options)`:
- it starts from that step's CLI defaults, read from its `_menu.py` parser (`_subcommand_defaults()`), so every `run()` parameter is present and the pipeline never drifts from `AAFTF <step>` defaults — do not hard-code step settings in `pipeline.py`
- `shared` pipeline options (`cpus`, `memory`, `workdir`, `debug`, `quiet`) go only to steps that have that option; a `None` value (e.g. no `--memory`) keeps the step default. Exception: `assemble` gets `{workdir}/assemble_{method}` (`_assemble_workdir()`), since SPAdes restarts from and MEGAHIT refuses an existing folder, and `filter` has already created the pipeline workdir
- `step_options` are the pipeline options specific to that step plus the file names chaining the steps; `pipe=True` is always set
- `tests/test_pipeline.py` checks, for all three pipelines, that each step receives exactly its CLI options, with defaults kept unless a pipeline option overrides them

## Pipeline Workflow

Recommended order for full assembly QC:

```
trim → [mito optional, PE only] → filter → assemble → vecscreen or fcs_screen
  → [taxpurge, needs a phylum] → covpurge → rmdup
  → [polish optional, mainly with long reads] → sort → assess → [depth optional]

separately, on the finished assemblies: [fcs_gx_purge, needs RAM for its database]
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
| sort | rmdup FASTA | `{base}.sorted.fasta` |
| depth | sorted FASTA + reads | `{base}.final.fasta` (MITO contigs tagged), `{base}.coverage_stats.txt`, `{base}.coverage_stats.mito_candidates.tsv` |
| assess | final FASTA | printed stats (optional report file) |

Run standalone, `polish` defaults to `<input prefix>.polished.fasta` and `depth` to `coverage_stats.txt`.

## External Tool Versions (minimum known-good)

- samtools >= 1.13 (pixi pins >= 1.24; needed for `sort --write-index` and `flagstat -O tsv`). No samtools version branching: aligner → sorted BAM pipelines use `utility.align_to_sorted_bam()`, which pipes SAM straight into `samtools sort` and writes the `.bai` index in the same step
- minimap2 >= 2.17
- mosdepth >= 0.3
- sourmash >= 4.8 with sourmash_plugin_branchwater >= 0.9.15 (`taxpurge`: `fastmultigather` + `tax genome`)
- SPAdes >= 3.15 for `assemble`
- NCBI FCS-adaptor 0.5.5 (hard-coded in `resources.py`)

## Code Style Guidelines

### Import Organization
- Standard library imports first, then third-party, then local AAFTF imports
- Use explicit imports (avoid `from module import *`)
- Group related imports together
- Import AAFTF modules with absolute imports (`from aaftf.utility import ...`)

```python
# Standard library
import os
import sys
import subprocess

# Third-party
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Local AAFTF imports
from aaftf.utility import check_file, run_cmd
from aaftf.resources import DB_LINKS
```

### Function and Variable Naming
- Use `snake_case` for function and variable names
- Use `UPPER_SNAKE_CASE` for constants
- Function names should be descriptive verbs (`run_assembly`, `check_file`, `calculate_stats`)
- Variable names should be clear and concise (`input_file`, `output_dir`, `cpu_count`)

### Documentation and Docstrings
- Use Google-style docstrings and type hints on all function signatures
- Every module should have a module docstring explaining its purpose
- Every public function should have a docstring with description, arguments, and returns

```python
def calculate_n50(contig_lengths):
    """Calculate N50 statistic from a list of contig lengths.

    Args:
        contig_lengths: List of integers representing contig lengths

    Returns:
        Integer N50 value
    """
```

### Error Handling
- Use try/except blocks for external tool execution
- Provide informative error messages with context
- Report progress with a module-level `logger = logging.getLogger(__name__)` (`logger.info`, `logger.warning`, `logger.error`); `aaftf.main.main()` configures it via `setup_logging()` from `-q/--quiet` and `-v/--verbose` (argparse dest `debug`). Don't add "ERROR:"/"WARNING:" prefixes or leading/trailing whitespace to messages — the formatter adds the level name
- Handle file existence checks with `check_file()` utility
- Never call `sys.exit()` in a subcommand or helper — raise instead, and let `aaftf.main.main()` log the message and return the exit code (the console script passes it to `sys.exit`): `FileNotFoundError` for a missing input file, tool or database (exit 2), `ValueError` for bad or missing arguments, `RuntimeError` when an external tool fails or produces no output (both exit 1); Ctrl-C returns 130. Put the whole message (including any hint) in the exception rather than logging it first, so it is not printed twice. The traceback is shown with `-v/--verbose`.

### Code Structure
- Package layout: `aaftf/` (one module per subcommand: `trim.py`, `assemble.py`, ...; `_menu.py` holds every subcommand's argparse parser; `main.py` is the entry point; `pipeline.py`; `utility.py`; `resources.py`; `data/`), `tests/` (pytest, `test_<module>.py`, shared fixtures in `conftest.py`), `docs/` (Sphinx)
- Common utilities are in `utility.py`
- Resource URLs and constants are in `resources.py`
- Each subcommand module has a `run(<CLI dest names>..., **kwargs)` function; its parser is a `<name>_menu()` in `aaftf/_menu.py` that binds it with `set_defaults(func=<module>.run)` (see "Architecture" above)

### CLI Argument Parsing
- Use `argparse` for command-line interface
- Follow existing patterns for common arguments (`-c/--cpus`, `-m/--memory`, `-o/--out`, `-1/--read1`, `-2/--read2`, `-w/--workdir`); path arguments use `type=str`
- Finish every subcommand's "optional arguments" group with `add_verbosity_args(optional)` (`-q/--quiet`, `-v/--verbose` with dest `debug`); `pipe` is a `run()`-only parameter, not a CLI option
- Provide help text for all arguments
- Use consistent naming across subcommands

### External Tool Integration
- Check tool availability with `require_tools([...])` (raises `FileNotFoundError` naming the missing tools); `AAFTF dependency` reports every tool
- Run external commands with `run_cmd()` / `execute()` (or `align_to_sorted_bam()` for aligner → sorted, indexed BAM)

### File I/O Patterns
- Support both compressed (.gz) and uncompressed files
- Use `check_file()` to validate input files
- Use `safe_remove()` for file cleanup. Create a subcommand's working directory with `workdir, custom_workdir = make_workdir(workdir, "<name>")` and remove it with `cleanup_workdir(workdir, debug, custom_workdir)` (it never deletes a user-supplied `--workdir`)
- Reference databases are listed in `aaftf.resources.DATABASES` and downloaded **only** by `AAFTF database`. Other subcommands get them with `require_databases([...])`, which returns the stored paths or raises `FileNotFoundError` naming the `AAFTF database ...` command to run — never download them from a subcommand. (Only `filter`'s `-a/--screen_accessions` / `-u/--screen_urls` fetch sequences themselves.)
- Every module follows the same top-level layout: (1) module docstring; (2) imports (stdlib, third-party, local; conditional `try: import ...` blocks right after them); (3) `__all__` listing the public constants, classes and functions; (4) constants and module-level variables, public first, then private ones and shared state (`logger`, `_cache`, ...); (5) public ABCs/Protocols; (6) public classes; (7) public functions — `run` first in subcommand modules, the rest in the order they are called; (8) private classes; (9) private helpers, also in call order; (10) the `if __name__ == "__main__":` block. Module-level code must not call anything defined further down (e.g. build objects lazily in a function, as `open_url()` does).
- Other shared helpers in `aaftf/utility.py` — use these instead of re-implementing: `download_file()` (never leaves a partial file behind; for `database` and filter's -a/-u), `db_file()` (database file location across `$AAFTF_DB` folders / `~/.cache/aaftf`), `require_tools()`, `run_cmd()` / `execute()`, `align_to_sorted_bam()`, `concat_files()`, `open_maybe_gz()`, `next_step_name()`, `basename_from_reads()`, `filter_fasta()` / `write_fasta()`, `calc_nx()`
- Handle file paths with `os.path` operations for cross-platform compatibility

### Logging and Output
- Use `logger.info/warning/error` for user-facing messages
- Use `print_cmd()` to show commands being executed
- Provide progress feedback for long-running operations
- Use debug flags for verbose output during development

### Performance Considerations
- Respect CPU and memory limits specified by user
- Use streaming operations for large files when possible
- Implement efficient file parsing with BioPython iterators
- Consider external tool memory requirements in parameter defaults

### Version Handling
- `aaftf/__init__.py` resolves `AAFTF.__version__` via `importlib.metadata.version("AAFTF")`, falling back to `"0.0.0+unknown"` if the package isn't installed
- The version string itself is set at build time by `hatch-vcs` from git tags (see `[tool.hatch.version]` in `pyproject.toml`)
- Main application (`main.py`) imports it with `from aaftf import __version__`
- Version is displayed both via `--version` flag and at application startup
- Maintains PEP 440 compatibility for Python packaging standards

### Dependencies and the bioconda recipe
- Dependencies are listed in several places that must stay in sync: `[project] dependencies` and the pixi features in `pyproject.toml` (then relock `pixi.lock`), `environment.yml` / `environment.dev.yml`, and the bioconda recipe `bioconda/meta.yaml`
- **Whenever a dependency is added, removed or its version bounds change (Python packages, external tools, or `requires-python`), update `bioconda/meta.yaml` in the same change**: its `run` requirements mirror `[project] dependencies` plus the default-pipeline tools (`[tool.pixi.dependencies]`) and `taxpurge`'s sourmash tools; optional tools are only named in `about.description`; `python_min` / `python_max` mirror `requires-python`
- bioconda's autobump bot only bumps the version and sha256 of the recipe in bioconda-recipes (`recipes/aaftf/meta.yaml`), never its dependencies, so dependency changes must also be pushed to bioconda-recipes by hand (e.g. onto the autobump PR for that release)
- Releases update the recipe automatically: for a stable release, `release.yml`'s `additional_bump_command` sets `version` and resets `build.number` to 0 in the release PR, and `publish.yml`'s `bioconda-sha256` job then opens a PR into `dev` with the PyPI sdist's `sha256`. Prereleases leave the recipe unchanged. For a recipe-only fix of the same version, bump `build.number` by hand instead

## Packaging, Containers and CI/Release

### Python package
- Built with hatchling + hatch-vcs (version from the `v*` git tag; no file carries it). The sdist is an allowlist (`[tool.hatch.build.targets.sdist] include`): `aaftf/`, `tests/`, `CHANGELOG.md`, the `patches/polca/` files, plus `pyproject.toml`, `README.md`, `LICENSE`. Add a path there only if building or testing the package needs it; docs, CI, container, dev files and `bioconda/` stay out
- The wheel installs `aaftf/` plus `aaftf_polca.sh` / `fix_consensus_from_vcf.pl` into `bin/` (hatch `shared-scripts`)

### Containers
- `Dockerfile` (pixi `complete` environment, source-built bowtie2/SPAdes/MEGAHIT for older CPUs) is the only published image: `ghcr.io/stajichlab/aaftf`. It writes `/etc/profile.d/aaftf.sh` so AAFTF stays on PATH in login shells (`bash -l`, used by SLURM/Nextflow) of a `.sif` converted from it; keep that hook when editing the Dockerfile
- No `.sif` is published: Singularity/Apptainer users run `apptainer pull aaftf.sif docker://ghcr.io/stajichlab/aaftf:<tag>`, which `docker-build.yml`'s `apptainer-test` job checks (including `bash -lc 'AAFTF --version'`). `AAFTF.def` stays for building locally (e.g. with databases built in); keep its `%post` steps in step with the Dockerfile
- Images are signed by the shared `container_build` (`sign: true`): only the index digest, stored as one `sha256-<index digest>` tag; buildx attaches the SLSA provenance. Verify with `cosign verify ghcr.io/stajichlab/aaftf:<tag> --certificate-identity-regexp '^https://github.com/stajichlab/gh_workflows/.github/workflows/container_build.yml@' --certificate-oidc-issuer https://token.actions.githubusercontent.com --certificate-github-workflow-repository stajichlab/AAFTF` (the certificate names the shared workflow, not AAFTF's)
- `container-cleanup.yml` (weekly): keeps only the latest alpha/beta tag, then the shared `container_cleanup` deletes untagged versions older than 14 days, keeping every manifest a kept image references, and removes orphaned signature tags. Never delete untagged GHCR versions with a plain age filter: the platform and attestation manifests of a tagged multi-arch image are untagged versions

### CI and release workflows
- Shared workflows come from `stajichlab/gh_workflows@v2` (`python_test_pip` / `python_test_pixi` / `python_test_conda`, `container_build`, `container_cleanup`, `release_prepare`, `release_publish`, `changelog_check`, `citation_check`, `sync_dev`, `pixi_lock_refresh`, `bioconda_recipe_solve`). A caller must grant every permission the shared job requests (e.g. `contents: read`), plus `id-token: write` for signing. There is no `setup_command`: put setup lines at the top of the multi-line `test_command`
- Releases: `release.yml` (manual; npm-style `bump_type`, e.g. `minor` turns `0.8.0-beta.0` into `0.8.0`, `preminor` starts `0.9.0-beta.0`) opens a `release/v<version>` PR from the `[Unreleased]` CHANGELOG section; merging it runs `publish.yml` (sdist build + tests, bioconda recipe solve, tag and GitHub release, PyPI, container, bioconda sha256 PR, sync into `dev`). Add a CHANGELOG `[Unreleased]` entry with every user-facing change
- There is no conda channel of our own: conda users install from bioconda (`recipe-solve.yml`, calling the shared `bioconda_recipe_solve`, on recipe/`pyproject.toml` changes and before each release, checks the recipe can be installed: a dry-run solve of the recipe's run requirements on every supported Python with glibc 2.17); see "Dependencies and the bioconda recipe"

## Testing Guidelines

Since AAFTF primarily integrates external bioinformatics tools, testing focuses on:
- pytest tests in `tests/`, marked `unit` (no external tools) or `integration` (needs the external tools)
- Mocked subprocess commands and files written (assert on side effects, not on the caller's kwargs)
- CLI argument parsing validation
- File format compatibility (FASTA/FASTQ parsing)

## Common Patterns

### Subcommand Module Structure
```python
"""Module description."""

import logging
from aaftf.utility import check_file, run_cmd

logger = logging.getLogger(__name__)


def run(input, debug=False, pipe=False, **kwargs):
    """Main entry point for subcommand."""
    # Validate inputs
    if not check_file(input):
        raise FileNotFoundError(f"Input file not found: {input}")

    logger.info("Starting processing...")

    # Execute external tools (prints the command; stderr shown only with -v/--verbose)
    run_cmd(["tool", "--input", input], debug)

    logger.info("Processing complete")
```

This project follows bioinformatics best practices with emphasis on reproducibility, proper tool integration, and clear user feedback.
