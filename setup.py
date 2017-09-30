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
