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

Most jobs call the org's shared reusable workflows from
[stajichlab/gh_workflows](https://github.com/stajichlab/gh_workflows) (pinned to `@v1`).

| Workflow | Runs on | What it does |
|---|---|---|
| `ci.yml` | every push to `main` and every pull request | pre-commit lint; unit tests on Python 3.10-3.14 (pip install, no external tools) and the full test suite, including integration tests, in the locked pixi `dev` environment (shared `python_test`); Sphinx docs build with warnings as errors |
| `ci-conda.yml` | pushes/PRs touching Python code or `environment.yml` | the full test suite in the `environment.yml` conda environment (shared `python_test`) |
| `changelog.yml` | every pull request | requires an entry under `[Unreleased]` in `CHANGELOG.md` (shared `changelog_check`; label `skip-changelog` to skip) |
| `citation.yml` | changes to `CITATION.cff` | validates it (shared `citation_check`) |
| `dependabot-pixi-lock.yml` | Dependabot PRs touching `pyproject.toml` | regenerates `pixi.lock` (shared `pixi_lock_refresh`) |
| `release.yml` | manual (`workflow_dispatch`) | computes the next version from the latest `v*` tag (npm-style: `major`/`minor`/`patch`, `pre*`, `prerelease`; e.g. `v0.7.0-beta.4` + `prerelease` -> `v0.7.0-beta.5`, + `patch` -> `v0.7.0`), moves the `[Unreleased]` notes in `CHANGELOG.md` to that version, updates `CITATION.cff` and opens a `Release v<version>` PR into `main` (shared `release_prepare`) |
| `publish.yml` | a merged `release/v*` PR | tags the release and creates the GitHub release (shared `release_publish`), then runs `python-publish.yml`, `conda-build.yml` and `docker-build.yml` for the tag and merges `main` back into `dev` (shared `sync_dev`) |
| `python-publish.yml` | called by `publish.yml` (or manual, with a tag) | builds the sdist/wheel from the tag, checks the built version matches it, uploads to PyPI and attaches them to the GitHub release |
| `conda-build.yml` | called by `publish.yml` (or manual) | builds the conda package from `ci/recipe/recipe.yaml` and uploads it to the `stajichlab` anaconda channel |
| `docker-build.yml` | called by `publish.yml`, or manual (`:main` test image) | builds and smoke-tests the Docker image (pixi `complete` environment) for ghcr.io (shared `container_build`), signs it with provenance, and attaches a Singularity `.sif` to releases |
| `container-cleanup.yml` | weekly | prunes untagged and old prerelease images from ghcr.io |

Add release notes under `## [Unreleased]` in `CHANGELOG.md` ([Keep a Changelog](https://keepachangelog.com) format) with every pull request. To release, run **Actions → Release**, review and merge the PR it opens. Set the optional `RELEASE_TOKEN` secret (a PAT or GitHub App token with contents and pull requests: write) so CI runs on the release PR; the repo must allow GitHub Actions to create pull requests (**Settings → Actions → General**). Run the CI checks locally with `pixi run -e dev pre-commit run --all-files` and `pixi run -e dev pytest tests/`.

## Code Style and Architecture

See [AGENTS.md](AGENTS.md) for full development guidelines, subcommand conventions, and code patterns.
