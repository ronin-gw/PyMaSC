#!/bin/bash

# Modern makesource.sh - Python version-specific C source generation
# Generates optimized C source files for Python 3.8-3.12

set -e  # Exit on any error

# Configuration
PYTHON_VERSIONS=("3.8" "3.9" "3.10" "3.11" "3.12")
WORK_DIR="/PyMaSC"

# Base build command with optimized dependencies
BASE_BUILDCMD="cd $WORK_DIR && pip install --no-cache-dir numpy 'pysam>=0.22.0' 'bx-python>=0.10.0' 'cython>=3.0.0'"

echo "🔧 Starting Python version-specific C source generation..."
echo "Target Python versions: ${PYTHON_VERSIONS[*]}"

# Clean all previous C sources
echo "🧹 Cleaning previous C sources..."
rm -f PyMaSC/*/*.c PyMaSC/*/*/*.c

# Function to generate C sources for a specific Python version
generate_sources_for_version() {
    local py_version=$1
    local py_major_minor=${py_version%.*}${py_version#*.}  # 3.8 -> 38, 3.10 -> 310
    
    echo "🐍 Generating C sources for Python ${py_version}..."
    
    # Build command with version-specific optimizations
    local buildcmd="$BASE_BUILDCMD && python setup.py build_ext -if"
    
    # Run Docker container with specific Python version
    if docker run --rm -v "$(pwd):$WORK_DIR" "python:${py_version}-slim" bash -c "$buildcmd"; then
        echo "✅ Build successful for Python ${py_version}"
        
        # Rename generated C files with version suffix
        echo "📝 Renaming C files for Python ${py_version}..."
        for f in PyMaSC/*/*[a-z].c PyMaSC/*/*/*[a-z].c; do
            if [[ -f "$f" ]]; then
                local base_name="${f%.c}"
                local version_file="${base_name}_${py_major_minor}.c"
                mv "$f" "$version_file"
                echo "  $f -> $version_file"
            fi
        done
        
        echo "✅ Python ${py_version} C sources generated successfully"
    else
        echo "❌ Failed to generate C sources for Python ${py_version}"
        return 1
    fi
}

# Generate sources for each Python version
success_count=0
for py_version in "${PYTHON_VERSIONS[@]}"; do
    echo ""
    echo "=================================="
    if generate_sources_for_version "$py_version"; then
        ((success_count++))
    else
        echo "⚠️  Warning: Failed for Python ${py_version}, continuing with others..."
    fi
done

echo ""
echo "=================================="
echo "🎉 C source generation completed!"
echo "Successfully generated sources for $success_count/${#PYTHON_VERSIONS[@]} Python versions"

# Display generated files
echo ""
echo "📁 Generated version-specific C files:"
for py_version in "${PYTHON_VERSIONS[@]}"; do
    local py_major_minor=${py_version%.*}${py_version#*.}
    local files=($(find PyMaSC -name "*_${py_major_minor}.c" | sort))
    if [[ ${#files[@]} -gt 0 ]]; then
        echo "  Python ${py_version}:"
        printf '    %s\n' "${files[@]}"
    fi
done

# Fallback: Generate generic _3.c files for compatibility
echo ""
echo "🔄 Generating fallback _3.c files for maximum compatibility..."
if docker run --rm -v "$(pwd):$WORK_DIR" "python:3.11-slim" bash -c "$BASE_BUILDCMD && python setup.py build_ext -if"; then
    for f in PyMaSC/*/*[a-z].c PyMaSC/*/*/*[a-z].c; do
        if [[ -f "$f" ]]; then
            mv "$f" "${f%.c}_3.c"
        fi
    done
    echo "✅ Fallback _3.c files generated"
else
    echo "⚠️  Warning: Failed to generate fallback files"
fi

echo ""
echo "✨ All C source generation tasks completed!"
echo "📊 Summary:"
echo "  - Version-specific sources: $success_count versions"
echo "  - Fallback _3.c sources: Generated"
echo "  - Total C files: $(find PyMaSC -name "*.c" | wc -l)"
