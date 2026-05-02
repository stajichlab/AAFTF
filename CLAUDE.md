# CLAUDE.md — Claude Code Project Context for AAFTF

See [AGENTS.md](AGENTS.md) for full development guidelines, code style, and common patterns.

## Key Architecture Reminders

- Entry point: `AAFTF/AAFTF_main.py` — defines all subcommand parsers and `run_subtool()` dispatcher
- Each subcommand is a module with a `run(parser, args)` function
- `AAFTF/pipeline.py` orchestrates the full end-to-end workflow using `Namespace` objects built from `args_dict`
- Aliases are registered via argparse `aliases=` and normalized by `ALIAS_MAP` in `run_subtool()`

## CLI Framework Conventions

- Every subcommand parser must have `-v/--debug` and `--pipe` flags — enforced by test suite (`TestAssessParser`, `TestSortParser`, `TestFixTblParser`)
- `ALIAS_MAP` in `AAFTF_main.py` maps all aliases to canonical names; update it when adding new aliases
- When building a `Namespace` for pipeline steps, include **all** attributes the target `run()` function accesses — check the submodule source to avoid `AttributeError`

## Subcommand Reference

| Canonical name | Aliases | Module | Key external tools |
|---|---|---|---|
| `trim` | `trim_reads`, `read_trim` | `AAFTF/trim.py` | bbduk.sh, trimmomatic, fastp |
| `filter` | `filter_reads`, `read_filter` | `AAFTF/filter.py` | bbduk.sh, bowtie2, bwa, minimap2, samtools |
| `assemble` | `asm`, `spades` | `AAFTF/assemble.py` | spades.py, megahit, unicycler |
| `vecscreen` | `vectorscreen`, `vector_blast` | `AAFTF/vecscreen.py` | blastn, makeblastdb |
| `fcs_screen` | `ncbi_fcs`, `ncbi_fcs-screen` | `AAFTF/fcs_screen.py` | run_fcsadaptor.sh (singularity/docker) |
| `fcs_gx_purge` | `ncbi_fcs-gx`, `ncbi_fcs_gx`, `gx` | `AAFTF/fcs_gx_purge.py` | run_gx.py (NCBI FCS-GX) |
| `sourpurge` | `purge` | `AAFTF/sourpurge.py` | sourmash, bwa, samtools |
| `rmdup` | `dedup` | `AAFTF/rmdup.py` | minimap2 |
| `polish` | `pilon`, `polca` | `AAFTF/polish.py` | pilon, nextPolish, polca, bwa, samtools |
| `sort` | — | `AAFTF/sort.py` | (none — BioPython only) |
| `assess` | `stats` | `AAFTF/assess.py` | (none — BioPython only) |
| `depth` | `coverage`, `cov` | `AAFTF/depth.py` | minimap2, bwa, samtools, mosdepth |
| `mito` | `mito_asm`, `mitochondria` | `AAFTF/mito.py` | NOVOPlasty, minimap2 |
| `fix_tbl` | `fix` | `AAFTF/fix_tbl.py` | (none) |
| `pipeline` | — | `AAFTF/pipeline.py` | all of the above |

## Module `args` Signatures

Key `args` attributes accessed by each `run()` function:

- **trim**: `left`, `right`, `basename`, `method`, `memory`, `cpus`, `minlen`, `avgqual`, `debug`, `pipe`, `trimmomatic`, `trimmomatic_adaptors`, `trimmomatic_leadingwindow`, `trimmomatic_trailingwindow`, `trimmomatic_slidingwindow`, `trimmomatic_quality`, `trimmomatic_clip`, `merge`, `dedup`, `cutfront`, `cuttail`, `cutright`
- **filter**: `workdir`, `cpus`, `AAFTF_DB`, `left`, `right`, `screen_accessions`, `screen_urls`, `screen_local`, `basename`, `aligner`, `memory`, `debug`, `pipe`
- **assemble**: `method`, `workdir`, `cpus`, `memory`, `isolate`, `careful`, `assembler_args`, `tmpdir`, `left`, `right`, `merged`, `pipe`, `debug`, `out`
- **vecscreen**: `workdir`, `AAFTF_DB`, `infile`, `outfile`, `cpus`, `percent_id`, `stringency`, `debug`, `pipe`
- **fcs_screen**: `AAFTF_DB`, `container_engine`, `workdir`, `infile`, `image`, `prok`, `fcs_script`, `outfile`, `debug`, `pipe`
- **fcs_gx_purge**: `workdir`, `db`, `input`, `taxid`, `outfile`, `debug`, `pipe`
- **sourpurge**: `workdir`, `cpus`, `left`, `right`, `sourdb`, `sourdb_type`, `input`, `AAFTF_DB`, `kmer`, `phylum`, `mincovpct`, `outfile`, `taxonomy`, `debug`, `pipe`
- **rmdup**: `workdir`, `cpus`, `input`, `percent_id`, `percent_cov`, `minlen`, `exhaustive`, `debug`, `out`, `pipe`
- **polish**: `method`, `memory`, `cpus`, `left`, `right`, `longreads`, `workdir`, `infile`, `outfile`, `iterations`, `debug`, `diploid`, `ploidy`, `pipe`, `polca`, `trimmomatic`
- **sort**: `input`, `minlen`, `out`, `name`
- **assess**: `input`, `report`, `telomere_monomer`, `telomere_n_repeat`, `telomere_window`
- **depth**: `input`, `out`, `left`, `right`, `longreads`, `illumina_preset`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `pipe`, `min_contig_len`, `no_plot`, `plot_format`, `quantize`, `quantize_labels`
- **mito**: `workdir`, `left`, `right`, `seed`, `reference`, `minlen`, `maxlen`, `out`, `starting`, `pipe`
- **fix_tbl**: `table`, `report`, `output`, `debug`, `pipe`

## Pipeline `create_namespace()` Convention

`pipeline.py` uses a `create_namespace(options, required_args=None, **extra)` helper:
- `options`: list of keys to copy from the pipeline `args_dict`
- `required_args`: dict of values to override/add (pipeline-step-specific settings)
- Always include `"pipe": True` in `required_args` for every step
- Always include `"debug"` in the `options` list so parent debug flag propagates
- Do NOT append attributes directly to the returned Namespace; include them in `required_args` before calling

## Pipeline Workflow

Recommended order for full assembly QC:

```
trim → [mito optional, PE only] → filter → assemble → vecscreen
  → [fcs_screen / fcs_gx_purge optional] → sourpurge → rmdup
  → polish → sort → assess → [depth optional]
```

Expected inputs/outputs per step:

| Step | Input | Output |
|---|---|---|
| trim | raw FASTQ | `{base}_1P.fastq.gz`, `{base}_2P.fastq.gz` |
| filter | trimmed FASTQ | `{base}_filtered_1.fastq.gz`, `{base}_filtered_2.fastq.gz` |
| assemble | filtered FASTQ | `{base}.{method}.fasta` |
| vecscreen | assembled FASTA | `{base}.vecscreen.fasta` |
| sourpurge | vecscreen FASTA | `{base}.sourpurge.fasta` |
| rmdup | sourpurge FASTA | `{base}.rmdup.fasta` |
| polish | rmdup FASTA | `{base}.polish.fasta` |
| sort | polished FASTA | `{base}.final.fasta` |
| assess | sorted FASTA | printed stats (optional report file) |
| depth | final FASTA + reads | `coverage_stats.txt` |

## External Tool Versions (minimum known-good)

- samtools >= 1.0 (many modules branch on samtools version)
- minimap2 >= 2.17
- mosdepth >= 0.3
- sourmash >= 4.x (LCA-based classification)
- SPAdes >= 3.15 for `assemble`
- NCBI FCS-adaptor 0.5.5 (hard-coded in `resources.py`)

## Recent Additions (2026-03-02)

### `depth` subtool — coverage analysis

**New file:** `AAFTF/depth.py`; registered in `AAFTF_main.py` with aliases `coverage`/`cov`.

**Workflow:**
1. Counts reads in each input FASTQ
2. Maps Illumina reads with `minimap2 -ax sr` (default) or `bwa mem` (`--aligner bwa`); maps long reads with `minimap2 -ax map-ont` (default) / `map-pb` / `map-hifi` (`--longread_type`)
3. Merges BAMs when both read types are provided (`samtools merge`)
4. Runs `samtools flagstat` per read type, then `mosdepth` on the combined BAM
5. Parses `mosdepth.summary.txt` for per-contig mean depths and `mosdepth.global.dist.txt` for coverage breadth
6. Flags contigs with mean depth > assembly_mean + 3×SD as possible contaminants/organelles (uses population SD; contigs are the full assembly population, not a sample)

**Key `args` attributes:** `input`, `out`, `left`, `right`, `longreads`, `illumina_preset`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `pipe`, `min_contig_len`, `no_plot`, `plot_format`, `quantize`, `quantize_labels`

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

### CLI consistency fixes (AAFTF_main.py)

- Added `--pipe` and `-v/--debug` to `fix_tbl` parser
- Added `--pipe` to `pipeline` parser
- Removed dead `runall` stub from `run_subtool()`
- Fixed assembler help text: "Default: unicycler" → "Default: spades"

### pipeline.py refactor

- `sort_args` and `assess_args` now use `create_namespace()` (consistent with all other steps)
- Assembly-specific attributes (`merged`, `isolate`, `careful`) now passed via `required_args` before `create_namespace()` call

### Bioinformatics enhancements

- `assess.py`: `findTelomere()` window is now configurable (`--telomere_window`, default 200 bp); added N-gap count and total N bases to stats output
- `rmdup.py`: now computes and reports both N50 and N75; filtering threshold uses N75 (as before)
- `depth.py`: added `--min_contig_len` (default 500 bp) to exclude short contigs from outlier depth analysis

## Running Tests

```bash
# Unit tests only (no external tools required)
conda run -n base python -m pytest tests/ -m unit -v

# All tests (integration tests require external tool environment)
conda run -n base python -m pytest tests/ -v

# Compile-check all modules
python -m py_compile AAFTF/filter.py AAFTF/polish.py AAFTF/trim.py \
    AAFTF/vecscreen.py AAFTF/sourpurge.py AAFTF/fcs_screen.py \
    AAFTF/depth.py AAFTF/AAFTF_main.py AAFTF/pipeline.py \
    AAFTF/assess.py AAFTF/rmdup.py
```
