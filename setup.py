#!/usr/bin/env python
import distutils
import subprocess
import os
from setuptools import setup, Extension
from setuptools.command import build_ext

BASEDIR = os.path.dirname(__file__)


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


def _setup():
    setup(
        cmdclass={
            "makeclean": MakeCleanCommand,
            "build_ext": BuildExtCommand
        },
        ext_modules=[
            Extension(
                "PyMaSC.reader.bigwig.bigwig",
                sources=[_basedir("PyMaSC/reader/bigwig/bigwig.pyx")],
                include_dirs=[_basedir("external/KentLib/inc")],
                extra_link_args=[_basedir("external/KentLib/lib/jkweb.a")]
            ),
            Extension(
                "PyMaSC.reader.bigwig.interval",
                sources=[_basedir("PyMaSC/reader/bigwig/interval.pyx")],
                include_dirs=[_basedir("external/KentLib/inc")],
                extra_link_args=[_basedir("external/KentLib/lib/jkweb.a")]
            )
        ],
        entry_points={
            "console_scripts": [
                "pymasc = PyMaSC.pymasc:exec_entrypoints",
                "calcmappablelen = PyMaSC.calcmappablelen:exec_entrypoints"
            ],
        }
    )


if __name__ == "__main__":
    _setup()
