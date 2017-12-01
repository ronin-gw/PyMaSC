#!/usr/bin/env python
import sys
import os
import shutil
from setuptools import setup, Extension, find_packages

from PyMaSC import VERSION

IMPORTERROR_FMT = "Failed to import {}. Install build time dependencies first " \
                  "(e.g. `pip install cython numpy pysam`) and try again."
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    sys.exit(IMPORTERROR_FMT.format("Cython"))
try:
    import numpy
except ImportError:
    sys.exit(IMPORTERROR_FMT.format("numpy"))
try:
    import pysam
except ImportError:
    sys.exit(IMPORTERROR_FMT.format("pysam"))

BASEDIR = os.path.dirname(__file__)

NUMPY_INCLUDE = numpy.get_include()
PYSAM_PATH = pysam.__path__[0]
PYSAM_INCLUDES = [PYSAM_PATH, os.path.join(PYSAM_PATH, "include", "htslib")]
sys.path.append(os.path.join(BASEDIR, "external", "bx-python", "lib"))

EXTRA_C_ARGS = ["-O3", "-ffast-math"]


def _basedir(path):
    return os.path.join(BASEDIR, path)


def _link_source():
    for ext in (".pyx", ".pxd"):
        shutil.copy2(_basedir("PyMaSC/utils/bx/binary_{}_endian".format(sys.byteorder) + ext),
                     _basedir("PyMaSC/utils/bx/binary_file" + ext))


def _define_extension(name, include_numpy=False):
    include_dirs = [NUMPY_INCLUDE] if include_numpy else []

    return Extension(
        name, sources=[_basedir(name.replace('.', '/') + ".pyx")],
        include_dirs=include_dirs, extra_compile_args=EXTRA_C_ARGS
    )


def _setup():
    #
    _link_source()

    #
    setup(
        name="PyMaSC",
        version=VERSION,
        author="Hayato Anzawa",
        author_email="anzawa@sb.ecei.tohoku.ac.jp",
        description="Python implementation to calc mappability-sensitive cross-correlation "
                    "for fragment length estimation and quality control for ChIP-Seq.",
        url="https://github.com/ronin-gw/PyMaSC",
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
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: Implementation :: CPython",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        keywords="NGS ChIP-Seq",
        platforms=["POSIX", "Mac OS X"],
        license="MIT",
        install_requires=[
            "numpy>=1.12.0",
            "pysam>=0.12.0.1",
            "scipy>=0.18.1",
            "bx-python>=0.7.3",
            "matplotlib>=2.0.0"
        ],
        packages=find_packages(),
        cmdclass={"build_ext": build_ext},
        ext_modules=cythonize([
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
        ] + [
            Extension(
                "PyMaSC.core.readlen",
                sources=[_basedir("PyMaSC/core/readlen.pyx")],
                include_dirs=PYSAM_INCLUDES,
                extra_compile_args=EXTRA_C_ARGS
            )
        ]),
        entry_points={
            "console_scripts": [
                "pymasc = PyMaSC.pymasc:exec_entrypoint",
                "pymasc-precalc = PyMaSC.calcmappablelen:exec_entrypoint"
            ],
        }
    )


if __name__ == "__main__":
    _setup()
