# PyMaSC CI/CD Workflows

This directory contains GitHub Actions workflows for continuous integration and testing of PyMaSC.

## Workflows

### 1. `build-check.yml` - Build Check
- **Trigger**: All pushes and PRs
- **Purpose**: Quick sanity check for build environment
- **Tests**: Dependencies, BitArray compilation, setup.py

### 2. `test.yml` - Basic Tests  
- **Trigger**: Push to master/modernize-python-support, PRs to master
- **Purpose**: Fast basic testing with Python 3.8 on Ubuntu
- **Tests**: Unit tests, imports, basic commands

### 3. `test-matrix.yml` - Matrix Tests
- **Trigger**: Push to master/modernize-python-support, PRs to master  
- **Purpose**: Test multiple Python versions (3.8-3.12) on Ubuntu and macOS
- **Tests**: Unit tests, quick integration tests, command-line tools

### 4. `test-full.yml` - Full Test Suite
- **Trigger**: Push to master, PRs to master, daily schedule, manual
- **Purpose**: Comprehensive testing including Golden tests
- **Tests**: All tests with coverage, numerical validation, Golden tests

## Test Strategy

1. **Build Check**: Runs on every push to catch basic issues early
2. **Basic Tests**: Quick feedback on core functionality  
3. **Matrix Tests**: Ensures compatibility across Python versions and OS
4. **Full Tests**: Complete validation including time-consuming tests

## Caching Strategy

- **pip cache**: Speeds up dependency installation
- **Test data cache**: Preserves large test files (BAM, BigWig)
- **Build artifacts**: Uploaded for debugging failed tests

## Known Issues

- Python 3.9+ requires using pre-compiled C sources due to Cython compatibility
- macOS requires x86_64 architecture flag for BitArray compilation
- Full test suite takes ~10-15 minutes due to Golden tests

## Local Testing

To run CI tests locally:

```bash
# Basic tests
pytest tests/unit/ -v

# Quick integration tests  
pytest tests/integration/ -k "not slow" -v

# Full test suite
pytest tests/ -v --cov=PyMaSC
```