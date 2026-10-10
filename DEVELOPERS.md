# Developer Guide

## Running Tests

Tests live in `tests/` and use pytest. Run from the repo root.

```bash
# Unit tests only — no external bioinformatics tools required
pixi run -e dev pytest tests/ -m unit -v

# All tests — integration tests require the full tool environment
pixi run -e dev pytest tests/ -v
```

## Compile Check

Quick syntax check of all modules without running anything:

```bash
python -m py_compile aaftf/*.py
```

## Lint

```bash
pixi run -e dev pre-commit run --all-files   # ruff, ruff-format, codespell, pydocstyle, ...
```

## Continuous Integration and Releases

GitHub Actions workflows in `.github/workflows/`:

| Workflow | Runs on | What it does |
|---|---|---|
| `ci.yml` | every push to `main` and every pull request | pre-commit lint; unit tests on Python 3.10-3.14 (pip install, no external tools); the full test suite, including integration tests, in the locked pixi `dev` environment; Sphinx docs build with warnings as errors |
| `publish_release.yml` | manual (`workflow_dispatch`) | bumps the version from the latest `v*` tag (npm-style: `major`/`minor`/`patch`, `pre*`, `prerelease`; e.g. `v0.7.0-beta.4` + `prerelease` -> `v0.7.0-beta.5`, + `patch` -> `v0.7.0`), updates `CITATION.cff`, turns the top `## ... (Development)` section of `CHANGES.md` into `## <version> (Stable)` or `(Pre-release)` and uses it as the release notes (stopping if it is missing or empty), adds a fresh `## Unreleased (Development)` section, runs the tests, tags, and creates the GitHub release |
| `python-publish.yml` | a published release | builds the sdist/wheel (full git history, so hatch-vcs picks up the tag), checks the built version matches the tag, uploads to PyPI |
| `conda-build.yml` | a published release (or manual) | builds the conda package from `ci/recipe/recipe.yaml` and uploads it to the `stajichlab` anaconda channel |
| `docker-build.yml` | pushes to `main`, published releases | builds the Docker image (pixi `complete` environment) for ghcr.io, signs it, and attaches a Singularity `.sif` to releases |
| `container-cleanup.yml` | weekly | prunes untagged and old prerelease images from ghcr.io |

Write release notes under the `## ... (Development)` heading at the top of `CHANGES.md` as you go. Run the CI checks locally with `pixi run -e dev pre-commit run --all-files` and `pixi run -e dev pytest tests/`.

## Code Style and Architecture

See [AGENTS.md](AGENTS.md) for full development guidelines, subcommand conventions, and code patterns.
