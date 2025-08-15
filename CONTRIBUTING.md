# Contributing and Release Guide

This document focuses on packaging and release flows. For development tasks (tests, type-checks, etc.), see the main `Makefile` and project README.

## Development Setup
- Python 3.8+ is required for development; 3.11+ is recommended for release tooling.
- Editable install with dev extras:
  - `make install-dev` or `pip install -e .[dev]`
- Build BitArray static library (needed for native extensions):
  - `make bitarray-lib` (clean with `make bitarray-clean`)

## Running Tests
- Full suite with coverage:
  - `make test` (HTML report under `htmlcov/`)
- Quick unit-only run:
  - `make test-quick`
- Specific groups:
  - `make test-unit`, `make test-integration`, `make test-all`
- Parallel execution (pytest-xdist):
  - `make test-parallel`, `make test-quick-parallel`, `make test-unit-parallel`, `make test-integration-parallel`, `make test-all-parallel`
  - Custom workers: `make test-workers WORKERS=4`
- Targeted runs:
  - Single file: `make test-<name>` (maps to `tests/unit/test_<name>.py`)
  - Pattern: `make test-match PATTERN='<expr>'`

## Type Checking
- Run mypy with repo config:
  - `make typecheck`

## Building Native Extensions (local)
- Build in place (incremental): `make build`
- Force rebuild everything: `make build-force`
- Full clean + rebuild: `make rebuild`
- Clean artifacts and caches: `make clean`

## Generating Versioned C Sources (Makefile.sources)
The project can pre-generate version-specific C sources (e.g., `*_310.c`) from Cython. This is mainly for sdist quality and cross-Python reproducibility; wheel builds via cibuildwheel do not require manual source generation.

Prerequisites
- Docker and Docker Compose are required for the containerized multi-Python build.

Common commands
- Recommended parallel build (isolated containers):
  - `make -f Makefile.sources parallel`
- Sequential build for all versions (3.8–3.13):
  - `make -f Makefile.sources all`
- Single-version builds:
  - `make -f Makefile.sources py38 | py39 | py310 | py311 | py312 | py313`
- Extract/sync C files from existing build containers:
  - `make -f Makefile.sources extract`
- Dependency/environment checks:
  - `make -f Makefile.sources check-deps` (includes BitArray/Cython/Docker checks)
  - `make -f Makefile.sources platform-info`
- Cleanup helpers:
  - Remove non-versioned C: `make -f Makefile.sources clean-non-versioned`
  - Remove versioned C: `make -f Makefile.sources clean-versioned`
  - Remove compiled `.so`: `make -f Makefile.sources clean-shared-libs`
  - Full cleanup (including Docker): `make -f Makefile.sources clean`

## Packaging Overview
- Wheels are built per-OS/architecture and per-CPython version due to C extensions (Cython/BitArray).
- We build:
  - Linux: manylinux wheels via cibuildwheel (Docker required)
  - macOS: universal2 (x86_64 + arm64) wheels, built natively on macOS
  - sdist: source distribution for fallback builds

## Prerequisites
- Python runner for cibuildwheel: Python 3.11+ (recommended)
- Docker installed (for Linux/manylinux wheels)
- macOS host for macOS wheels

We keep the cibuildwheel toolchain isolated from dev dependencies via an extra named `release`.

```
# Create a Python 3.11+ virtualenv (example with pyenv)
pyenv install -s 3.11.9
pyenv virtualenv 3.11.9 pymasc-release
pyenv activate pymasc-release

# Install minimal release toolchain
pip install -e .[release]
```

## Build Commands (Release)
Use the dedicated `Makefile.release` so the release flow is independent from the dev Makefile.

- Source distribution:
  - `make -f Makefile.release sdist`
- Linux wheels (manylinux):
  - `make -f Makefile.release wheels-linux`
- macOS wheels (universal2):
  - `make -f Makefile.release wheels-macos`
- Both:
  - `make -f Makefile.release wheels-all`
- List artifacts:
  - `make -f Makefile.release check`

Notes:
- You can override the runner Python if not active, e.g.:
  - `make -f Makefile.release RUNPY=python3.11 wheels-all`
- Build outputs are placed under `dist/`.

## Upload to PyPI
1. Ensure `~/.pypirc` is set or use `-u/-p` options interactively.
2. Run:
   - `python -m twine upload dist/*`

Tip: use TestPyPI first:
```
python -m twine upload --repository testpypi dist/*
```

## Matrix and Configuration
- Python versions: 3.8–3.13
- Linux arches: manylinux x86_64 by default (aarch64 optional via `CIBW_ARCHS_LINUX=aarch64`)
- macOS: universal2 (x86_64 + arm64)
- Configuration lives in `pyproject.toml` under `[tool.cibuildwheel]`.

Override examples:
```
# Only cp311 and cp312
CIBW_BUILD="cp311-* cp312-*" make -f Makefile.release wheels-linux
# Add aarch64 for Linux (QEMU; slower)
CIBW_ARCHS_LINUX=aarch64 make -f Makefile.release wheels-linux
```

### Linux arch selection (CIBW_ARCHS_LINUX)
- Default: builds `x86_64` only.
- x86_64 only: `CIBW_ARCHS_LINUX=x86_64 make -f Makefile.release wheels-linux`
- aarch64 only: `CIBW_ARCHS_LINUX=aarch64 make -f Makefile.release wheels-linux`
- Both at once: `CIBW_ARCHS_LINUX="x86_64 aarch64" make -f Makefile.release wheels-linux`

Notes
- aarch64 builds use QEMU emulation and are slower.
- Docker is required (uses manylinux images).
- Per-command overrides are fine; no need to persist environment variables.

## C Sources vs Wheels
- `Makefile.sources` + `docker-compose.build.yml` are for generating versioned C files from Cython.
- For wheels, cibuildwheel installs build-time deps (Cython, numpy, pysam, pyBigWig) as configured in `pyproject.toml` and compiles the extensions during wheel build.
- You generally do not need to run `Makefile.sources` to build wheels.

## Troubleshooting
- Runner Python warning: cibuildwheel v3 requires Python 3.11+. Use a 3.11+ runner (this does not limit wheel target versions).
- BitArray build flags:
  - Linux: handled by `[tool.cibuildwheel.linux].before-build`
  - macOS: universal2 flags are set in `[tool.cibuildwheel.macos].before-build`
- If building fails due to environment issues, try a clean workspace: `make clean`.
