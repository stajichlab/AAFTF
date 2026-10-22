[![CI](https://github.com/stajichlab/AAFTF/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/stajichlab/AAFTF/actions/workflows/ci.yml)
[![CI/conda](https://github.com/stajichlab/AAFTF/actions/workflows/ci-conda.yml/badge.svg?branch=main)](https://github.com/stajichlab/AAFTF/actions/workflows/ci-conda.yml)
[![Python](https://img.shields.io/badge/python-3.10_%7C_3.11_%7C_3.12_%7C_3.13_%7C_3.14-blue?logo=python)](https://github.com/stajichlab/AAFTF/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/stajichlab/AAFTF/graph/badge.svg)](https://codecov.io/gh/stajichlab/AAFTF)
[![License](https://img.shields.io/github/license/stajichlab/AAFTF?label=license)](https://github.com/stajichlab/AAFTF/blob/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/AAFTF?label=PyPI)](https://pypi.org/project/AAFTF/)
[![Container](https://ghcr-badge.egpl.dev/stajichlab/aaftf/latest_tag?trim=major&label=Image)](https://github.com/stajichlab/AAFTF/pkgs/container/aaftf)
[![Documentation Status](https://app.readthedocs.org/projects/aaftf/badge/?version=latest)](https://aaftf.readthedocs.io/en/latest/)

# AAFTF - Automatic Assembly For The Fungi
*Authors: Jason Stajich, Jon Palmer and Cheng-Hung Tsai*

![AAFTF logo](docs/AAFTF.png)

AAFTF automates draft genome assembly from Illumina short reads (optionally with long reads), with
read trimming and contaminant filtering, assembly, vector and contaminant screening of contigs,
duplicate removal, contig sorting and renaming, and summary statistics including telomere detection.
Each step is a subcommand that can be run on its own, or all together with `AAFTF pipeline_short` (Illumina), `pipeline_long` (ONT/PacBio) or `pipeline_hybrid`.

**Documentation:** <https://aaftf.readthedocs.io/en/latest/>

# Install
With conda (bioconda) or [pixi](https://pixi.sh), from a checkout of this repository:

```
$ conda env create -f environment.yml && conda activate aaftf   # what the pipelines need
$ pixi install -e complete && pixi shell -e complete             # or every optional tool too
$ AAFTF dependency                                               # check which tools are found
```

Docker and Singularity images, the pixi environments and the list of external tools are described in
[Installation](https://aaftf.readthedocs.io/en/latest/installation.html).

Then download the reference databases once (into `~/.cache/aaftf`, or the folder(s) in `$AAFTF_DB`):

```
$ AAFTF database required     # what the default pipelines need
$ AAFTF database taxpurge     # sourmash databases for taxpurge (pipelines with -p)
```

# Quick start

```
# Illumina reads
AAFTF pipeline_short -1 STRAINX_R1.fq.gz -2 STRAINX_R2.fq.gz -o STRAINX -p Ascomycota -c 16 -m 64

# long reads only
AAFTF pipeline_long -lr STRAINX_ont.fq.gz -o STRAINX --genome_size 40m -c 16

# Illumina + long reads
AAFTF pipeline_hybrid -1 STRAINX_R1.fq.gz -2 STRAINX_R2.fq.gz -lr STRAINX_ont.fq.gz -o STRAINX -c 16 -m 64
```

Each writes `STRAINX.final.fasta` and prints the assembly statistics. `-p/--phylum` is optional:
without it the taxonomic screen (`taxpurge`) is skipped. Steps whose output already exists are skipped,
so an interrupted run can simply be restarted. See
[Workflow](https://aaftf.readthedocs.io/en/latest/workflow.html) for each pipeline's steps and files,
and for running the steps one by one.

# Commands
`AAFTF -h` lists the subcommands; `AAFTF <command> -h` shows its options. Each is documented in
[Commands](https://aaftf.readthedocs.io/en/latest/commands/index.html).

| Group | Subcommands |
|---|---|
| Setup | `dependency`, `database` |
| Steps | `trim`, `mito`, `filter`, `assemble`, `vecscreen`, `taxpurge`, `fcs_screen`, `covpurge`, `rmdup`, `polish`, `sort`, `assess`, `depth` |
| Pipeline | `pipeline_short`, `pipeline_long`, `pipeline_hybrid` |
| Utilities | `fcs_gx_purge`, `fix_tbl` |

# Notes
This is partially a Python rewrite of [JAAWS](https://github.com/nextgenusfs/jaaws), a Unix shell based
cleanup and assembly tool written by Jon. Changes between versions are listed in [CHANGES.md](CHANGES.md).

# Authors
* Jason Stajich [@hyphaltip](https://github.com/hyphaltip) - http://lab.stajich.org, [@hyphaltip.bsky](https://bsky.app/profile/hyphaltip.bsky.social)
* Jon Palmer [@nextgenusfs](https://github.com/nextgenusfs) - [@jonpalmer.bsky](https://bsky.app/profile/jonpalmer.bsky.social)
* Cheng-Hung Tsai

# Citation
Palmer JM and Stajich JE. (2023). Automatic assembly for the fungi (AAFTF): genome assembly pipeline (v0.5.0). Zenodo. doi: 10.5281/zenodo.1620526
