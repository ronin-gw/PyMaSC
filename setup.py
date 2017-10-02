#!/usr/bin/env python
import sys
import os
import shutil
from setuptools import setup, Extension

import numpy
from Cython.Build import cythonize

BASEDIR = os.path.dirname(__file__)
NUMPY_INCLUDE = numpy.get_include()
EXTRA_C_ARGS = ["-O3", "-ffast-math"]

sys.path.append(os.path.join(BASEDIR, "external", "bx-python", "lib"))


def _basedir(path):
    return os.path.join(BASEDIR, path)


def _define_extension(name, include_numpy=False):
    include_dirs = [NUMPY_INCLUDE] if include_numpy else []

    return Extension(
        name, sources=[_basedir(name.replace('.', '/') + ".pyx")],
        include_dirs=include_dirs, extra_compile_args=EXTRA_C_ARGS
    )


def _setup():
    #
    for ext in ".pyx", ".pxd":
        shutil.copy2(_basedir("PyMaSC/utils/bx/binary_{}_endian".format(sys.byteorder) + ext),
                     _basedir("PyMaSC/utils/bx/binary_file" + ext))

    #
    setup(
        name="PyMaSC",
        version="0.1.0",
        author="Hayato Anzawa",
        author_email="anzawa@sb.ecei.tohoku.ac.jp",
        description="Python implementation to calc mappability-sensitive cross-correlation "
                    "for fragment length estimation and quality control for ChIP-Seq.",
        download_url="https://github.com/ronin-gw/PyMaSC",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: POSIX",
            "Programming Language :: Cython",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: Implementation :: CPython",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        platforms=["POSIX", "Mac OS X"],
        license="MIT",
        install_requires=[
            "numpy>=1.12.0",
            "pysam>=0.12.0.1",
            "scipy>=0.18.1",
            "bx-python>=0.7.3",
            "matplotlib>=2.0.0"
        ],
        ext_modules=cythonize(
            [_define_extension(name, include_numpy=True) for name in (
                "PyMaSC.reader.bx.bbi_file",
                "PyMaSC.reader.bx.bigwig_file",
                "PyMaSC.reader.bigwig",
                "PyMaSC.core.mappability",
                "PyMaSC.core.ncc",
                "PyMaSC.core.mscc"
            )] + [_define_extension(name) for name in (
                "PyMaSC.reader.bx.bpt_file",
                "PyMaSC.core.readlen",
                "PyMaSC.utils.bx.binary_file"
            )]
        ),
        entry_points={
            "console_scripts": [
                "pymasc = PyMaSC.pymasc:exec_entrypoint",
                "pymasc-precalc = PyMaSC.calcmappablelen:exec_entrypoint"
            ],
        }
    )


if __name__ == "__main__":
    _setup()
