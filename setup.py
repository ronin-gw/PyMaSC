#!/usr/bin/env python
import sys
import os
import shutil
from setuptools import setup, Extension

import numpy
import pysam

BASEDIR = os.path.dirname(__file__)
NUMPY_INCLUDE = numpy.get_include()
PYSAM_PATH = pysam.__path__[0]
PYSAM_INCLUDES = [PYSAM_PATH, os.path.join(PYSAM_PATH, "include", "htslib")]
EXTRA_C_ARGS = ["-O3", "-ffast-math"]

try:
    CYTHON_BUILD = True
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
    sys.path.append(os.path.join(BASEDIR, "external", "bx-python", "lib"))
except ImportError:
    CYTHON_BUILD = False


def _basedir(path):
    return os.path.join(BASEDIR, path)


def _link_source():
    exts = (".pyx", ".pxd") if CYTHON_BUILD else (".c", )

    for ext in exts:
        shutil.copy2(_basedir("PyMaSC/utils/bx/binary_{}_endian".format(sys.byteorder) + ext),
                     _basedir("PyMaSC/utils/bx/binary_file" + ext))


def _define_extension(name, include_numpy=False):
    ext = ".pyx" if CYTHON_BUILD else ".c"
    include_dirs = [NUMPY_INCLUDE] if include_numpy else []

    return Extension(
        name, sources=[_basedir(name.replace('.', '/') + ext)],
        include_dirs=include_dirs, extra_compile_args=EXTRA_C_ARGS
    )


def _get_setup_args():
    cmdclass = {} if not CYTHON_BUILD else {"build_ext": build_ext}
    ext_modules = [
        _define_extension(name, include_numpy=True) for name in (
            "PyMaSC.reader.bx.bbi_file",
            "PyMaSC.reader.bx.bigwig_file",
            "PyMaSC.reader.bigwig",
            "PyMaSC.core.mappability",
            "PyMaSC.core.ncc",
            "PyMaSC.core.mscc"
        )
    ] + [
        _define_extension(name) for name in (
            "PyMaSC.reader.bx.bpt_file",
            "PyMaSC.utils.bx.binary_file"
        )
    ] + [Extension(
        "PyMaSC.core.readlen",
        sources=[_basedir("PyMaSC/core/readlen" + (".pyx" if CYTHON_BUILD else ".c"))],
        include_dirs=PYSAM_INCLUDES,
        extra_compile_args=EXTRA_C_ARGS
    )]

    return cmdclass, ext_modules


def _setup():
    #
    _link_source()
    cmdclass, ext_modules = _get_setup_args()
    if CYTHON_BUILD:
        ext_modules = cythonize(ext_modules)

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
        cmdclass=cmdclass,
        ext_modules=ext_modules,
        entry_points={
            "console_scripts": [
                "pymasc = PyMaSC.pymasc:exec_entrypoint",
                "pymasc-precalc = PyMaSC.calcmappablelen:exec_entrypoint"
            ],
        }
    )


if __name__ == "__main__":
    _setup()
