#!/usr/bin/env python
import distutils
import subprocess
import os
import sys
import shutil
from setuptools import setup, Extension
from setuptools.command import build_ext

import numpy
from Cython.Build import cythonize

BASEDIR = os.path.dirname(__file__)
NUMPY_INCLUDE = numpy.get_include()
EXTRA_C_ARGS = ["-O3", "-ffast-math"]

sys.path.append(os.path.join(BASEDIR, "external", "bx-python", "lib"))


def _basedir(path):
    return os.path.join(BASEDIR, path)


class MakeCleanCommand(distutils.cmd.Command):
    description = "Run `make clean` for external libraries."
    user_options = []

    def initialize_options(self):
        pass

    def run(self):
        command = ' '.join(["cd", _basedir("external/KentLib"), "&&", "make", "clean"])
        self.announce("Clean KentLib ...", level=distutils.log.INFO)
        subprocess.check_call(command, shell=True)

    def finalize_options(self):
        pass


class BuildExtCommand(build_ext.build_ext):
    def run(self):
        self.run_command("makeclean")
        make_command = ' '.join(["cd", _basedir("external/KentLib"), "&&", "make"])
        subprocess.check_call(make_command, shell=True)
        build_ext.build_ext.run(self)


def _make_link():
    for ext in ".pyx", ".pxd":
        shutil.copy2(_basedir("PyMaSC/utils/bx/binary_{}_endian".format(sys.byteorder) + ext),
                     _basedir("PyMaSC/utils/bx/binary_file" + ext))


def _define_extension(name, include_numpy=False):
    include_dirs = [NUMPY_INCLUDE] if include_numpy else []

    return Extension(
        name, sources=[_basedir(name.replace('.', '/') + ".pyx")],
        include_dirs=include_dirs, extra_compile_args=EXTRA_C_ARGS
    )


def _setup():
    setup(
        cmdclass={
            "makeclean": MakeCleanCommand,
            "build_ext": BuildExtCommand
        },
        ext_modules=[_define_extension(name, include_numpy=True) for name in (
            "PyMaSC.reader.bx.bbi_file",
            "PyMaSC.reader.bx.bigwig_file",
            "PyMaSC.reader.bigwig",
            "PyMaSC.core.mappability",
            "PyMaSC.core.ncc",
            "PyMaSC.core.mscc"
        )] + cythonize([
            "PyMaSC/reader/bx/bpt_file.pyx",
            "PyMaSC/core/readlen.pyx",
            "PyMaSC/utils/bx/binary_file.pyx"
        ]),
        entry_points={
            "console_scripts": [
                "pymasc = PyMaSC.pymasc:exec_entrypoint",
                "pymasc-precalc = PyMaSC.calcmappablelen:exec_entrypoint"
            ],
        }
    )


if __name__ == "__main__":
    _make_link()
    _setup()
