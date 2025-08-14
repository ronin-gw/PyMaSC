#!/usr/bin/env python
"""
PyMaSC Setup Script
Modern setuptools-based setup with pyproject.toml integration
"""

import sys
import os
import subprocess
import glob
import platform
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command import build_ext

# Import build-time dependencies
import numpy
import pysam
import logging

# Configure logging to show source selection details
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

#
BASEDIR = Path(__file__).parent.absolute()

# Advanced source selection with staged fallback
# Priority: .pyx (Cython) -> version-specific C -> error


def _get_source_file(base_name):
    """
    Advanced source file selection with staged fallback logic.

    Environment Variables:
    - CYTHON_VERSION_SUFFIX: Direct versioned filename generation (e.g., "_38")
    - BUILD_MODE: "versioned" enables direct versioned file generation

    Fallback order:
    1. .pyx files (if Cython is available) - for development builds
    2. Version-specific C sources (*_39.c) - for optimized distribution
    3. Error if none found
    """
    py_version = f"{sys.version_info.major}{sys.version_info.minor}"

    # Check for environment variable-controlled versioned filename generation
    version_suffix = os.environ.get('CYTHON_VERSION_SUFFIX')
    build_mode = os.environ.get('BUILD_MODE')

    candidates = []

    # If in versioned build mode and we have a version suffix, generate versioned filenames directly
    if build_mode == "versioned" and version_suffix:
        logging.info(f"Versioned build mode detected: suffix={version_suffix}")
        # For versioned builds, check if we already have the versioned C file
        versioned_c_file = f"{base_name}{version_suffix}.c"
        if Path(versioned_c_file).exists():
            # Use existing versioned C file
            candidates.append((versioned_c_file, "Existing versioned C source"))
        else:
            # Need to build from .pyx first, then use versioned C file
            if CYTHON_AVAILABLE:
                candidates.append((f"{base_name}.pyx", "Cython source for versioned generation"))
            # Also add the target versioned C file in case it gets generated during build
            candidates.append((versioned_c_file, "Target versioned C source"))
    else:
        # Standard fallback logic
        candidates.extend([
            (f"{base_name}.pyx", "Cython source"),
            (f"{base_name}_{py_version}.c", f"Python {sys.version_info.major}.{sys.version_info.minor} optimized C source")
        ])

    for candidate_file, description in candidates:
        if Path(candidate_file).exists():
            logging.info(f"Using {description}: {candidate_file}")
            return candidate_file

    # No suitable source found
    available_files = []
    for pattern in [f"{base_name}.pyx", f"{base_name}_*.c"]:
        available_files.extend(glob.glob(pattern))

    raise FileNotFoundError(
        f"No suitable source file found for {base_name}\n"
        f"Environment: BUILD_MODE={build_mode}, CYTHON_VERSION_SUFFIX={version_suffix}\n"
        f"Searched for: {[f for f, _ in candidates]}\n"
        f"Available files: {available_files}"
    )


# Determine if Cython is available
CYTHON_AVAILABLE = False
try:
    from Cython.Build import cythonize
    CYTHON_AVAILABLE = True
    logging.info("Cython is available - .pyx files will be prioritized")
except ImportError:
    cythonize = lambda x: x  # noqa: ignore
    logging.info("Cython not available - using pre-built C sources")
# numpy path
NUMPY_INCLUDES = [numpy.get_include()]
# pysam path
PYSAM_PATH = pysam.__path__[0]
PYSAM_INCLUDES = [PYSAM_PATH, str(Path(PYSAM_PATH) / "include" / "htslib")]
# bx-python path removed - using pyBigWig implementation
# BitArray path
BITARRAY_DIR = str(BASEDIR / "external" / "BitArray")
BITARRAY_INCLUDES = [
    BITARRAY_DIR,
    str(BASEDIR / "external")
]

EXTRA_C_ARGS = ["-O3", "-ffast-math"]
# ARM64 Mac compatible flags
EXTRA_BA_ARGS = ["-Wextra", "-Wc++-compat"]

# Platform-specific linker flags for BitArray
BITARRAY_LIB_PATH = str(Path("external") / "BitArray" / "libbitarr.a")
BITARRAY_OBJ_PATH = str(Path("external") / "BitArray" / "bit_array.o")
if platform.system() == "Darwin":
    # macOS uses -all_load
    BITARRAY_LINK_ARGS = ["-Wl,-all_load", BITARRAY_LIB_PATH]
else:
    # Linux: link object file directly to avoid archive issues
    BITARRAY_LINK_ARGS = [BITARRAY_OBJ_PATH]


# _basedir function removed - using direct path handling in _get_source_file


class BuildExtCommand(build_ext.build_ext):
    def run(self):
        # Use relative path for BitArray directory
        bitarray_rel_path = str(Path("external") / "BitArray")
        command = ' '.join(["cd", bitarray_rel_path, "&&", "make", "libbitarr.a"])
        subprocess.check_call(command, shell=True)

        # Handle versioned filename generation for Cython builds
        version_suffix = os.environ.get('CYTHON_VERSION_SUFFIX')
        build_mode = os.environ.get('BUILD_MODE')

        if build_mode == "versioned" and version_suffix and CYTHON_AVAILABLE:
            logging.info(f"Versioned build mode: will generate files with suffix {version_suffix}")
            self._generate_versioned_sources(version_suffix)
            # Update extension sources to use versioned files
            self._update_extensions_for_versioned_sources(version_suffix)

        build_ext.build_ext.run(self)

    def _generate_versioned_sources(self, version_suffix):
        """
        Generate versioned C source files from .pyx files.
        This handles the custom versioned filename generation.
        """

        # Define the modules that need versioned C files
        core_modules = [
            "PyMaSC/reader/bx/bbi_file",
            "PyMaSC/reader/bx/bigwig_file",
            "PyMaSC/reader/bigwig",
            "PyMaSC/core/mappability",
            "PyMaSC/core/successive/ncc",
            "PyMaSC/core/successive/mscc",
            "PyMaSC/core/bitarray/bitarray",
            "PyMaSC/core/bitarray/mscc",
            "PyMaSC/core/readlen"
        ]

        for module_path in core_modules:
            pyx_file = f"{module_path}.pyx"
            target_c_file = f"{module_path}{version_suffix}.c"

            if Path(pyx_file).exists():
                logging.info(f"Generating versioned C source: {pyx_file} -> {target_c_file}")

                try:
                    # Generate C source with versioned filename directly
                    self._cythonize_to_versioned_file(pyx_file, target_c_file)
                    logging.info(f"  Created: {target_c_file}")

                except Exception as e:
                    logging.error(f"  Failed to generate {target_c_file}: {e}")
                    continue

    def _cythonize_to_versioned_file(self, pyx_file, target_c_file):
        """
        Generate versioned C file with proper race condition protection.
        Uses unique temporary filename within the source directory to prevent conflicts
        while preserving directory structure for Cython imports.
        """
        import shutil
        from Cython.Build import cythonize

        # Generate unique temporary C file name using process ID
        base_path = pyx_file[:-4]  # Remove .pyx extension
        temp_c_file = f"{base_path}_temp_{os.getpid()}.c"

        try:
            # Use Cython with build_dir=None to build in source directory
            # but temporarily rename the expected output to avoid conflicts
            _ = cythonize([pyx_file],
                          compiler_directives={'language_level': 3},
                          quiet=True,
                          build_dir=None)

            # Standard expected C file location
            standard_c_file = pyx_file.replace('.pyx', '.c')

            if Path(standard_c_file).exists():
                # Remove target file if it exists
                if Path(target_c_file).exists():
                    os.remove(target_c_file)
                # Move from standard location to versioned name
                shutil.move(standard_c_file, target_c_file)
            else:
                raise FileNotFoundError(f"Cython did not generate expected file: {standard_c_file}")

        finally:
            # Clean up any temporary files that might have been created
            for temp_pattern in [temp_c_file, f"{base_path}.c"]:
                if Path(temp_pattern).exists() and temp_pattern != target_c_file:
                    try:
                        os.remove(temp_pattern)
                    except OSError:
                        pass  # File might already be gone

    def _update_extensions_for_versioned_sources(self, version_suffix):
        """
        Update extension module sources to use versioned C files.
        This is called after versioned C files are generated.
        """
        logging.info(f"Updating extension sources to use versioned files with suffix {version_suffix}")

        # Update each extension to use the versioned C file
        updated_count = 0
        for ext in self.extensions:
            if len(ext.sources) == 1:
                source_file = ext.sources[0]
                # Check if this is a .pyx file that should be replaced with versioned C
                if source_file.endswith('.pyx'):
                    base_name = source_file[:-4]  # Remove .pyx extension
                    versioned_c_file = f"{base_name}{version_suffix}.c"

                    if Path(versioned_c_file).exists():
                        ext.sources[0] = versioned_c_file
                        logging.info(f"  Updated {ext.name}: {source_file} -> {versioned_c_file}")
                        updated_count += 1
                    else:
                        logging.warning(f"  Versioned C file not found: {versioned_c_file}")
                elif source_file.endswith('.c') and not source_file.endswith(f'{version_suffix}.c'):
                    # Check if there's a versioned version of this C file
                    base_name = source_file[:-2]  # Remove .c extension
                    versioned_c_file = f"{base_name}{version_suffix}.c"

                    if Path(versioned_c_file).exists():
                        ext.sources[0] = versioned_c_file
                        logging.info(f"  Updated {ext.name}: {source_file} -> {versioned_c_file}")
                        updated_count += 1

        logging.info(f"Updated {updated_count} extensions to use versioned C sources")


def _define_extension(name, include_dirs=[], extra_link_args=[], extra_compile_args=[]):
    """
    Define an extension module with advanced source file selection.
    Uses staged fallback: .pyx -> version-specific C -> generic C
    """
    base_path = name.replace('.', '/')
    source_file = _get_source_file(base_path)

    return Extension(
        name,
        sources=[source_file],
        include_dirs=include_dirs,
        extra_link_args=extra_link_args,
        extra_compile_args=EXTRA_C_ARGS + extra_compile_args
    )


def _build_extensions():
    """
    Build extension modules with proper dependency handling.
    Uses Cython compilation only if .pyx sources are being used.
    """
    extensions = [
        # numpy dependent modules
        _define_extension(name, include_dirs=NUMPY_INCLUDES) for name in (
            "PyMaSC.reader.bigwig",
            "PyMaSC.core.mappability",
            "PyMaSC.core.successive.ncc",
            "PyMaSC.core.successive.mscc"
        )
    ] + [
        # BitArray dependent modules
        _define_extension(
            name,
            include_dirs=BITARRAY_INCLUDES + NUMPY_INCLUDES,
            extra_link_args=BITARRAY_LINK_ARGS,
            extra_compile_args=EXTRA_BA_ARGS
        ) for name in (
            "PyMaSC.core.bitarray.bitarray",
            "PyMaSC.core.bitarray.mscc"
        )
    ] + [
        # pysam dependent modules
        _define_extension(name, include_dirs=PYSAM_INCLUDES) for name in (
            "PyMaSC.core.readlen",
        )
    ]

    # Apply Cython compilation only if .pyx sources are being used
    if CYTHON_AVAILABLE:
        # Check if any extension is using .pyx sources
        has_pyx_sources = any(
            any(source.endswith('.pyx') for source in ext.sources)
            for ext in extensions
        )

        if has_pyx_sources:
            logging.info("Compiling .pyx sources with Cython")
            return cythonize(extensions)
        else:
            logging.info("Using pre-built C sources (no .pyx files found)")
    else:
        logging.info("Using pre-built C sources (Cython not available)")

    return extensions


def _setup():
    """
    Main setup function - configuration moved to pyproject.toml
    """
    setup(
        cmdclass={"build_ext": BuildExtCommand},
        ext_modules=_build_extensions()
    )


if __name__ == "__main__":
    _setup()
