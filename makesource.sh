#!/bin/bash

# Modern makesource.sh - Docker Compose based C source generation
# Generates version-specific C source files for Python 3.8-3.12

set -e  # Exit on any error

# Configuration
PYTHON_VERSIONS=("3.8" "3.9" "3.10" "3.11" "3.12")

echo "🔧 Starting Python version-specific C source generation..."
echo "Target Python versions: ${PYTHON_VERSIONS[*]}"

# Function to rename generated C sources for a specific Python version
rename_sources_for_version() {
    local py_version=$1
    local py_major_minor=${py_version%.*}${py_version#*.}  # 3.8 -> 38, 3.10 -> 310

    echo "📝 Renaming C files for Python ${py_version}..."

    # Find and rename all generated C files
    for f in PyMaSC/*/*[a-z].c PyMaSC/*/*/*[a-z].c; do
        if [[ -f "$f" ]]; then
            local base_name="${f%.c}"
            local version_file="${base_name}_${py_major_minor}.c"
            mv "$f" "$version_file"
            echo "  $f -> $version_file"
        fi
    done
}

# Generate C sources for each Python version
success_count=0

for py_version in "${PYTHON_VERSIONS[@]}"; do
    echo ""
    echo "=================================="
    echo "🐍 Processing Python ${py_version}..."

    # Clean only non-versioned C sources (preserve *_XX.c files)
    echo "🧹 Cleaning non-versioned C sources..."
    find PyMaSC -name "*.c" ! -name "*_[0-9][0-9].c" ! -name "*_[0-9][0-9][0-9].c" -delete

    # Run specific service
    service_name="build-py${py_version//./}"  # 3.8 -> build-py38

    echo "🔨 Building C sources for Python ${py_version}..."
    if docker-compose -f docker-compose.build.yml up --build "$service_name"; then
        echo "✅ Build successful for Python ${py_version}"

        # Rename generated C files with version suffix
        rename_sources_for_version "$py_version"

        ((success_count++))
    else
        echo "❌ Failed to generate C sources for Python ${py_version}"
        echo "⚠️  Warning: Failed for Python ${py_version}, continuing with others..."
    fi

    # Clean up after each version
    echo "🧹 Cleaning up Docker containers..."
    docker-compose -f docker-compose.build.yml down -v
done

# Generate fallback _3.c files for maximum compatibility
echo ""
echo "🔄 Generating fallback _3.c files for maximum compatibility..."

# Clean and generate generic files using Python 3.11 as base
rm -f PyMaSC/*/*.c PyMaSC/*/*/*.c

if docker-compose -f docker-compose.build.yml up --build build-py311; then
    for f in PyMaSC/*/*[a-z].c PyMaSC/*/*/*[a-z].c; do
        if [[ -f "$f" ]]; then
            mv "$f" "${f%.c}_3.c"
        fi
    done
    echo "✅ Fallback _3.c files generated"
else
    echo "⚠️  Warning: Failed to generate fallback files"
fi

# Final cleanup
docker-compose -f docker-compose.build.yml down -v

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

# Display fallback files
echo ""
echo "📁 Generated fallback _3.c files:"
local fallback_files=($(find PyMaSC -name "*_3.c" | sort))
if [[ ${#fallback_files[@]} -gt 0 ]]; then
    printf '  %s\n' "${fallback_files[@]}"
fi

echo ""
echo "✨ All C source generation tasks completed!"
echo "📊 Summary:"
echo "  - Version-specific sources: $success_count versions"
echo "  - Fallback _3.c sources: Generated"
echo "  - Total C files: $(find PyMaSC -name "*.c" | wc -l)"