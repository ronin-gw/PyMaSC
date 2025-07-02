#!/bin/bash
set -e

# Script to test PyMaSC with multiple Python versions
PYTHON_VERSIONS=("3.8.9" "3.9.21" "3.10.15" "3.11.1" "3.12.11" "3.13.3")
FAILED_VERSIONS=()

echo "=== Testing PyMaSC with multiple Python versions ==="
echo

for version in "${PYTHON_VERSIONS[@]}"; do
    echo "================================"
    echo "Testing with Python $version"
    echo "================================"
    
    # Skip if version not installed
    if ! pyenv versions | grep -q "  $version"; then
        echo "Python $version not installed, skipping..."
        continue
    fi
    
    # Set python version
    export PYENV_VERSION=$version
    
    # Show python version
    python --version
    
    # Clean previous builds
    echo "Cleaning previous builds..."
    rm -rf build/ *.egg-info
    find . -name "*.so" -delete
    find . -name "*.o" -delete
    
    # Install dependencies
    echo "Installing dependencies..."
    pip install --quiet --upgrade pip
    pip install --quiet numpy "pysam>=0.22.0" "bx-python>=0.10.0"
    pip install --quiet cython pytest
    
    # Build BitArray
    echo "Building BitArray..."
    cd external/BitArray
    make clean >/dev/null 2>&1 || true
    make CC="gcc -arch x86_64" libbitarr.a
    cd ../..
    
    # Build PyMaSC
    echo "Building PyMaSC..."
    python setup.py build_ext --inplace
    
    # Run basic tests
    echo "Running basic tests..."
    if python -m pytest tests/unit/test_imports.py -v; then
        echo "✅ Python $version: Basic tests passed"
        
        # Run more comprehensive tests for stable versions
        if [[ "$version" == "3.8.9" || "$version" == "3.9.21" || "$version" == "3.10.15" ]]; then
            echo "Running comprehensive tests..."
            if python -m pytest tests/unit/ -v --tb=short; then
                echo "✅ Python $version: Unit tests passed"
            else
                echo "❌ Python $version: Unit tests failed"
                FAILED_VERSIONS+=("$version")
            fi
        fi
    else
        echo "❌ Python $version: Basic tests failed"
        FAILED_VERSIONS+=("$version")
    fi
    
    echo
done

# Reset to default Python version
pyenv shell 3.8.9

# Summary
echo "================================"
echo "Summary"
echo "================================"
echo "Tested Python versions: ${PYTHON_VERSIONS[@]}"
if [ ${#FAILED_VERSIONS[@]} -eq 0 ]; then
    echo "✅ All versions passed!"
else
    echo "❌ Failed versions: ${FAILED_VERSIONS[@]}"
    exit 1
fi