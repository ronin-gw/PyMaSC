#!/usr/bin/env python
import sys
import os
import subprocess
from setuptools import setup, Extension, find_packages
from setuptools.command import build_ext

from PyMaSC import VERSION


# Check build-time dependencies
NOTFOUNDMODULES = []
try:
    import numpy
except ImportError:
    NOTFOUNDMODULES.append("numpy")
try:
    import pysam
except ImportError:
    NOTFOUNDMODULES.append("pysam==0.15.1")

if NOTFOUNDMODULES:
    msg = "Failed to import build time dependencies.\nPlease install {} first (e.g. `pip install {}`) and try again.".format(
        " and ".join(NOTFOUNDMODULES),
        ' '.join(NOTFOUNDMODULES)
    )

    try:
        import Cython  # noqa
        msg += "\ncython is found in your environtment. If you do not want to build from cython source, remove cython before build."
    except ImportError:
        msg += "\nIf you want to build from cython source, install cython before build additionally."

    sys.exit(msg)


#
BASEDIR = os.path.abspath(os.path.dirname(__file__))

# cython
EXTENSION_SUFFIX = ".pyx"
try:
    from Cython.Build import cythonize
except ImportError:
    EXTENSION_SUFFIX = "_{}.c".format(sys.version_info[0])
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
EXTRA_BA_ARGS = ["-Wextra", "-Wc++-compat", "-march=native"]


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


def _setup():
    #
    setup(
        name="PyMaSC",
        version=VERSION,
        author="Hayato Anzawa",
        author_email="anzawa@sb.ecei.tohoku.ac.jp",
        description="Python implementation to calc mappability-sensitive cross-correlation "
                    "for fragment length estimation and quality control for ChIP-Seq.",
        long_description='\n'.join((
            "To estimate mean fragment length for single-ended sequencing data, "
            "cross-correlation between positive- and negative-strand reads is "
            "commonly used. One of the problems with this approach is phantom "
            "peak, which is the occasionally observed peak corresponding to the "
            "read length. In the ChIP-Seq guidelines by ENCODE consortia, cross-"
            "correlation at fragment length and read length are used for quality "
            "control metrics. Additionally, library length itself is one of the "
            "important parameters for analysises. However, estimating correct "
            "flagment length is not a easy task because of phantom peak occurrence.",
            "P Ramachandran et al. proposed MaSC, mappability-sensitive cross-"
            "correlation to remove the bias caused from ununiformity of "
            "mappability throughout the genome. This method provides cross-"
            "correlation landscape without phantom peak and much accurate mean "
            "fragment length estimation.",
            "PyMaSC is a tool implemented by python and cython to visualize "
            "(mappability-sensitive) cross-correlation and estimate ChIP-Seq "
            "quality metrics and mean fragment length with MaSC algorithm."
        )),
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
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: Implementation :: CPython",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        keywords="NGS ChIP-Seq",
        platforms=["POSIX", "Mac OS X"],
        license="MIT",
        install_requires=[
            "numpy>=1.12.0",
            "pysam==0.15.1",
            "bx-python>=0.7.3",
            "matplotlib>=2.0.0"
        ],
        packages=find_packages(),
        cmdclass={"build_ext": BuildExtCommand},
        ext_modules=cythonize([
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
                extra_link_args=[os.path.join(BITARRAY_DIR, "libbitarr.a")],
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
        ]),
        entry_points={
            "console_scripts": [
                "pymasc = PyMaSC.pymasc:main",
                "pymasc-precalc = PyMaSC.calcmappablelen:main",
                "pymasc-plot = PyMaSC.plot:main"
            ],
        }
    )


if __name__ == "__main__":
    _setup()
