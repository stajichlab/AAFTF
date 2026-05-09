# Developer Guide

## Running Tests

Tests live in `tests/` and use pytest. Run from the repo root.

```bash
# Unit tests only — no external bioinformatics tools required
python -m pytest tests/ -m unit -v

# All tests — integration tests require the full tool environment
python -m pytest tests/ -v
```

## Compile Check

Quick syntax check of all modules without running anything:

```bash
python -m py_compile AAFTF/filter.py AAFTF/polish.py AAFTF/trim.py \
    AAFTF/vecscreen.py AAFTF/sourpurge.py AAFTF/fcs_screen.py \
    AAFTF/depth.py AAFTF/AAFTF_main.py AAFTF/pipeline.py \
    AAFTF/assess.py AAFTF/rmdup.py AAFTF/download.py
```

## Code Style and Architecture

See [AGENTS.md](AGENTS.md) for full development guidelines, subcommand conventions, and code patterns.
