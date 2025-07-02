#!/usr/bin/env python
"""
PyMaSC Setup Script
Modern setuptools-based setup with pyproject.toml integration
"""

import sys
import os
import subprocess
from setuptools import setup, Extension
from setuptools.command import build_ext

# Import build-time dependencies
try:
    import numpy
except ImportError:
    raise ImportError("numpy is required for building PyMaSC. Please install it first: pip install numpy")

try:
    import pysam
except ImportError:
    raise ImportError("pysam is required for building PyMaSC. Please install it first: pip install 'pysam>=0.22.0'")

#
BASEDIR = os.path.abspath(os.path.dirname(__file__))

# cython
EXTENSION_SUFFIX = ".pyx"
try:
    from Cython.Build import cythonize
except ImportError:
    EXTENSION_SUFFIX = "_3.c"
    cythonize = lambda x: x
# numpy path
NUMPY_INCLUDES = [numpy.get_include()]
# pysam path
PYSAM_PATH = pysam.__path__[0]
PYSAM_INCLUDES = [PYSAM_PATH, os.path.join(PYSAM_PATH, "include", "htslib")]
# bx-python path
sys.path.append(os.path.join(BASEDIR, "external", "bx-python", "lib"))
# BitArray path
BITARRAY_DIR = os.path.join(BASEDIR, "external", "BitArray")
BITARRAY_INCLUDES = [
    BITARRAY_DIR,
    os.path.join(BASEDIR, "external")
]

EXTRA_C_ARGS = ["-O3", "-ffast-math"]
# ARM64 Mac compatible flags
EXTRA_BA_ARGS = ["-Wextra", "-Wc++-compat"]


def _basedir(path):
    return os.path.join(BASEDIR, path)


class BuildExtCommand(build_ext.build_ext):
    def run(self):
        command = ' '.join(["cd", BITARRAY_DIR, "&&", "make", "libbitarr.a"])
        subprocess.check_call(command, shell=True)
        build_ext.build_ext.run(self)


def _define_extension(name, include_dirs=[], extra_link_args=[], extra_compile_args=[]):
    return Extension(
        name,
        sources=[_basedir(name.replace('.', '/') + EXTENSION_SUFFIX)],
        include_dirs=include_dirs,
        extra_link_args=extra_link_args,
        extra_compile_args=EXTRA_C_ARGS + extra_compile_args
    )


def _build_extensions():
    """
    Build extension modules with proper dependency handling
    """
    return cythonize([
        # numpy dependent modules
        _define_extension(name, include_dirs=NUMPY_INCLUDES) for name in (
            "PyMaSC.reader.bx.bbi_file",
            "PyMaSC.reader.bx.bigwig_file", 
            "PyMaSC.reader.bigwig",
            "PyMaSC.core.mappability",
            "PyMaSC.core.ncc",
            "PyMaSC.core.mscc"
        )
    ] + [
        # BitArray dependent modules
        _define_extension(
            name,
            include_dirs=BITARRAY_INCLUDES,
            extra_link_args=["-Wl,-all_load", os.path.join(BITARRAY_DIR, "libbitarr.a")],
            extra_compile_args=EXTRA_BA_ARGS
        ) for name in (
            "PyMaSC.bacore.bitarray",
            "PyMaSC.bacore.mscc"
        )
    ] + [
        # pysam dependent modules
        _define_extension(name, include_dirs=PYSAM_INCLUDES) for name in (
            "PyMaSC.core.readlen",
        )
    ] + [
        # others
        _define_extension(name) for name in (
            "PyMaSC.reader.bx.bpt_file",
            "PyMaSC.utils.bx.binary_file"
        )
    ])


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