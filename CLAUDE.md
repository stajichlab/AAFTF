# CLAUDE.md — Claude Code Project Context for AAFTF

See [AGENTS.md](AGENTS.md) for full development guidelines, code style, and common patterns.

## Key Architecture Reminders

- Entry point: `AAFTF/AAFTF_main.py` — defines all subcommand parsers and `run_subtool()` dispatcher
- Each subcommand is a module with a `run(parser, args)` function
- `AAFTF/pipeline.py` orchestrates the full end-to-end workflow using `Namespace` objects built from `args_dict`
- Aliases are registered via argparse `aliases=` and normalized by `ALIAS_MAP` in `run_subtool()`

## CLI Framework Conventions

- Every subcommand parser should have `-v/--debug` and `--pipe` flags
- `ALIAS_MAP` in `AAFTF_main.py` maps all aliases to canonical names; update it when adding new aliases
- When building a `Namespace` for pipeline steps, include **all** attributes the target `run()` function accesses — check the submodule source to avoid `AttributeError`

## Recent Additions (2026-03-02)

### `depth` subtool — coverage analysis

**New file:** `AAFTF/depth.py`; registered in `AAFTF_main.py` with aliases `coverage`/`cov`.

**Workflow:**
1. Counts reads in each input FASTQ
2. Maps Illumina reads with `minimap2 -ax sr` (default) or `bwa mem` (`--aligner bwa`); maps long reads with `minimap2 -ax map-ont` (default) / `map-pb` / `map-hifi` (`--longread_type`)
3. Merges BAMs when both read types are provided (`samtools merge`)
4. Runs `samtools flagstat` per read type, then `mosdepth` on the combined BAM
5. Parses `mosdepth.summary.txt` for per-contig mean depths and `mosdepth.global.dist.txt` for coverage breadth
6. Flags contigs with mean depth > assembly_mean + 3×SD as possible contaminants/organelles

**Key `args` attributes accessed by `depth.run()`:**
`input`, `out`, `left`, `right`, `longreads`, `illumina_preset`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `pipe`

## Recent Fixes (2026-03-02)

### Prompt used to generate plan

> Investigate and fix CLI framework inconsistencies in AAFTF: runtime bugs in `pipeline.py`, broken alias routing in `run_subtool()`, and missing `--pipe`/`--debug` flags on some subcommand parsers.

### Changes implemented

**`AAFTF/pipeline.py`**
- Fixed `assess_args` Namespace to include `telomere_monomer` and `telomere_n_repeat` (prevented `AttributeError` at runtime)
- Removed `asm_args.spades_tmpdir = None` — wrong attribute name; `tmpdir` already extracted via `assembleOpts`

**`AAFTF/AAFTF_main.py`**
- Added `ALIAS_MAP` dict and alias normalization in `run_subtool()` so all argparse aliases route correctly
- Consolidated redundant `elif` branches (`pilon`, `fix`) into the map
- Added `-v/--debug` to `parser_mito`, `parser_sort`, `parser_assess`
- Added `--pipe` to `parser_sort`, `parser_assess`
- Added `-v` shorthand to `parser_rmdup`'s existing `--debug` flag
