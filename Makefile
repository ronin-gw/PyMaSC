# PyMaSC Development Makefile
# This is a convenience wrapper for common development tasks
# For C source generation from Cython files, see Makefile.sources

.PHONY: test test-quick test-all test-integration test-unit build install install-dev install-test install-plot install-requirements clean sources sources-parallel help

# Default target
help:
	@echo "PyMaSC Development Commands:"
	@echo "  make test          - Run full test suite with coverage"
	@echo "  make test-quick    - Run unit tests only (fast)"
	@echo "  make test-unit     - Run unit tests verbosely"
	@echo "  make test-integration - Run integration tests"
	@echo "  make test-all      - Run all tests verbosely"
	@echo "  make build         - Build Cython extensions in-place"
	@echo "  make install       - Install package in development mode"
	@echo "  make install-dev   - Install with all development dependencies"
	@echo "  make install-test  - Install with testing dependencies"
	@echo "  make install-plot  - Install with plotting dependencies"
	@echo "  make install-requirements - Install minimal runtime dependencies"
	@echo "  make clean         - Remove build artifacts and caches"
	@echo "  make sources       - Generate C sources from .pyx files"
	@echo "  make sources-parallel - Generate C sources in parallel"
	@echo ""
	@echo "For detailed C source generation options:"
	@echo "  make -f Makefile.sources help"

# Run tests with coverage (default test command)
test:
	python -m pytest tests/ --verbose --tb=short --cov=PyMaSC --cov-report=term-missing --cov-report=html

# Quick test run (unit tests only, no coverage)
test-quick:
	python -m pytest tests/unit/ -v

# Run unit tests verbosely
test-unit:
	python -m pytest tests/unit/ -v

# Run integration tests
test-integration:
	python -m pytest tests/integration/ -v

# Run all tests verbosely
test-all:
	python -m pytest tests/ -v

# Build Cython extensions in-place
build:
	python setup.py build_ext --inplace

# Install package in development mode
install:
	pip install -e .

# Install for development with all dependencies
install-dev:
	pip install -e .[dev]

# Install for testing (without plot dependencies)
install-test:
	pip install -e .[test]

# Install for plotting features
install-plot:
	pip install -e .[plot]

# Install using requirements.txt (minimal runtime dependencies)
install-requirements:
	pip install -r requirements.txt

# Clean build artifacts and caches
clean:
	# Remove Python caches
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	find . -type f -name "*.pyo" -delete 2>/dev/null || true
	
	# Remove build directories
	rm -rf build/ dist/ *.egg-info
	
	# Remove coverage data
	rm -rf htmlcov/ .coverage .coverage.*
	
	# Remove pytest cache
	rm -rf .pytest_cache/
	
	# Remove compiled Cython extensions
	find PyMaSC -name "*.so" -delete 2>/dev/null || true
	# Remove C files but preserve version-specific ones (*_38.c, *_39.c, *_310.c, *_311.c, *_312.c, *_3.c)
	find PyMaSC -name "*.c" ! -name "*_38.c" ! -name "*_39.c" ! -name "*_310.c" ! -name "*_311.c" ! -name "*_312.c" ! -name "*_3.c" -delete 2>/dev/null || true
	
	@echo "Clean complete."

# Generate C sources from Cython files
sources:
	$(MAKE) -f Makefile.sources all

# Generate C sources in parallel
sources-parallel:
	$(MAKE) -f Makefile.sources parallel

# Additional useful targets

# Run specific test file
test-%:
	python -m pytest tests/unit/test_$*.py -v

# Run tests matching a pattern
test-match:
	@if [ -z "$(PATTERN)" ]; then \
		echo "Usage: make test-match PATTERN=<pattern>"; \
		echo "Example: make test-match PATTERN=golden"; \
		exit 1; \
	fi
	python -m pytest tests/ -k "$(PATTERN)" -v

# Check code style (if linting is configured)
lint:
	@echo "Note: Linting not configured. Add flake8/black/isort as needed."

# Run type checking (if configured)
typecheck:
	@echo "Note: Type checking not configured. Add mypy as needed."